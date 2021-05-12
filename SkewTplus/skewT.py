"""
Module based on the skewT.py module from the SkewT Python project developed by
Thomas Chubb and the Matplotlib's examples/specialty_plots/skewt.py example.
 
In this version, all the SkewT plotting capabilities are now part of the 
SkewX axis class. This make the SkewT axis compatible with 
Matplolib's pyplot framework.


A new figure class with it a Figure Manager incorporated for easy visualization
and saving of plots is included.
The behavior of the pyplot **show** method is now part of the Figure class.
"""

# For python 3 portability
from __future__ import absolute_import, division, print_function, unicode_literals

import matplotlib.axis as maxis
import matplotlib.spines as mspines
import matplotlib.transforms as transforms
import numpy
from contextlib import ExitStack
from matplotlib import pyplot
from matplotlib.axes import Axes
from matplotlib.axes import Axes
from matplotlib.projections import register_projection
from matplotlib.projections import register_projection
from matplotlib.ticker import MultipleLocator, ScalarFormatter
from numpy import (
    array,
    linspace,
    log,
    where,
    concatenate,
    log10,
    logspace,
    zeros,
    logical_or,
)
from numpy.ma.core import masked_invalid, masked_array, getmaskarray

from SkewTplus.thermodynamics import (
    degCtoK,
    Rs_da,
    Cp_da,
    moistAscent,
    liftParcel,
    parcelAnalysis,
    virtualTemp3,
    virtualTemp4,
)


# The sole purpose of this class is to look at the upper, lower, or total
# interval as appropriate and see what parts of the tick to draw, if any.
class SkewXTick(maxis.XTick):
    def draw(self, renderer):
        # When adding the callbacks with `stack.callback`, we fetch the current
        # visibility state of the artist with `get_visible`; the ExitStack will
        # restore these states (`set_visible`) at the end of the block (after
        # the draw).
        with ExitStack() as stack:
            for artist in [
                self.gridline,
                self.tick1line,
                self.tick2line,
                self.label1,
                self.label2,
            ]:
                stack.callback(artist.set_visible, artist.get_visible())
            needs_lower = transforms.interval_contains(
                self.axes.lower_xlim, self.get_loc()
            )
            needs_upper = transforms.interval_contains(
                self.axes.upper_xlim, self.get_loc()
            )
            self.tick1line.set_visible(self.tick1line.get_visible() and needs_lower)
            self.label1.set_visible(self.label1.get_visible() and needs_lower)
            self.tick2line.set_visible(self.tick2line.get_visible() and needs_upper)
            self.label2.set_visible(self.label2.get_visible() and needs_upper)
            super().draw(renderer)

    def get_view_interval(self):
        return self.axes.xaxis.get_view_interval()


# This class exists to provide two separate sets of intervals to the tick,
# as well as create instances of the custom tick
class SkewXAxis(maxis.XAxis):
    def _get_tick(self, major):
        return SkewXTick(self.axes, None, major=major)

    def get_view_interval(self):
        return self.axes.upper_xlim[0], self.axes.lower_xlim[1]


# This class exists to calculate the separate data range of the
# upper X-axis and draw the spine there. It also provides this range
# to the X-axis artist for ticking and gridlines
class SkewSpine(mspines.Spine):
    def _adjust_location(self):
        pts = self._path.vertices
        if self.spine_type == "top":
            pts[:, 0] = self.axes.upper_xlim
        else:
            pts[:, 0] = self.axes.lower_xlim


# This class handles registration of the skew-xaxes as a projection as well
# as setting up the appropriate transformations. It also overrides standard
# spines and axes instances as appropriate.
class SkewXAxes(Axes):
    """
    Class to handle Skew-X axis used in the SkewT plots.

    This class contain all the necessary functions to plot vertical
    soundings with its corresponding parcel analysis.

    This class derives from the Matplolib's Axes_ class.

    For axis initialization the following keywords are supported to control the
    SkewT diagram properties:


    .. _Axes: https://matplotlib.org/api/axes_api.html#matplotlib.axes.Axes


    Parameters
    ----------
    pmin : float, optional
        Top level pressure value in hPa
    pmax : float, optional
        Bottom level pressure value in hPa
    tmax : float, optional
        skew X axis maximum temperature in Celsius
    tmin : float, optional
    skew X axis minimum temperature in Celsius

    Other Parameters
    ----------------

    kwargs : extra keywords arguments , optional
        See Axes_ class initialization parameters for more information on the the supported keywords
        The xscale and yscale has keywords have no effect.


    """

    # The projection must specify a name.  This will be used be the
    # user to select the projection, i.e.
    # ``subplot(111, projection='skewx')``.
    name = "skewx"

    def _init_axis(self):
        # Taken from Axes and modified to use our modified X-axis
        self.xaxis = SkewXAxis(self)
        self.spines["top"].register_axis(self.xaxis)
        self.spines["bottom"].register_axis(self.xaxis)
        self.yaxis = maxis.YAxis(self)
        self.spines["left"].register_axis(self.yaxis)
        self.spines["right"].register_axis(self.yaxis)

    def _gen_axes_spines(self):
        spines = {
            "top": SkewSpine.linear_spine(self, "top"),
            "bottom": mspines.Spine.linear_spine(self, "bottom"),
            "left": mspines.Spine.linear_spine(self, "left"),
            "right": mspines.Spine.linear_spine(self, "right"),
        }
        return spines

    def _set_lim_and_transforms(self):
        """
        This is called once when the plot is created to set up all the
        transforms for the data, text and grids.
        """
        rot = 30

        # Get the standard transform setup from the Axes base class
        super()._set_lim_and_transforms()

        # Need to put the skew in the middle, after the scale and limits,
        # but before the transAxes. This way, the skew is done in Axes
        # coordinates thus performing the transform around the proper origin
        # We keep the pre-transAxes transform around for other users, like the
        # spines for finding bounds
        self.transDataToAxes = (
                self.transScale + self.transLimits + transforms.Affine2D().skew_deg(rot,
                                                                                    0)
        )
        # Create the full transform from Data to Pixels
        self.transData = self.transDataToAxes + self.transAxes

        # Blended transforms like this need to have the skewing applied using
        # both axes, in axes coords like before.
        self._xaxis_transform = (
                transforms.blended_transform_factory(
                    self.transScale + self.transLimits, transforms.IdentityTransform()
                )
                + transforms.Affine2D().skew_deg(rot, 0)
                + self.transAxes
        )

    @property
    def lower_xlim(self):
        return self.axes.viewLim.intervalx

    @property
    def upper_xlim(self):
        pts = [[0.0, 1.0], [1.0, 1.0]]
        return self.transDataToAxes.inverted().transform(pts)[:, 0]

    def __init__(self, *args, **kwargs):
        """ New constructor """

        pmin = kwargs.pop("pmin", 100)
        pmax = kwargs.pop("pmax", 1050)
        tmax = kwargs.pop("tmax", 50)
        tmin = kwargs.pop("tmin", -30)
        _ = kwargs.pop("yscale", None)
        _ = kwargs.pop("xscale", None)

        super().__init__(yscale="log", xscale="linear", *args, **kwargs)

        self.setLimits(tmin=tmin, tmax=tmax, pmin=pmin, pmax=pmax)

        self.other_housekeeping()

        P = logspace(log10(pmax), log10(pmin), 100)

        self.linesHandlers = list()
        self.linesLabels = list()

        w = array(
            [
                0.00001,
                0.0001,
                0.0004,
                0.001,
                0.002,
                0.004,
                0.007,
                0.01,
                0.016,
                0.024,
                0.032,
            ]
        )
        self.add_mixratio_isopleths(
            w, P[P >= 700], color="g", ls="--", alpha=1.0, lw=0.5
        )

        self.add_dry_adiabats(
            linspace(250, 500, 18) - degCtoK, P, color="g", ls="--", alpha=1.0, lw=0.5
        )

        self.add_moist_adiabats(
            linspace(0, 44, 12), pmax, color="g", ls="--", alpha=1.0, lw=0.5
        )

        self.set_title("Sounding")

    def setLimits(self, tmin=-30, tmax=50, pmin=100, pmax=1000):
        """
        Set the limits of the Skew-T diagram

        Parameters
        ----------

        tmin : float
            Minimum temperature in celsius (skew x-Axis)

        tmax : float
            Maximum temperature in celsius (skew x-Axis)

        pmin : float
            Minimum pressure in hPa (log y-Axis)

        pmax : float
            Maximum pressure in hPa (log y-Axis)
        """

        self.pmin = pmin
        self.pmax = pmax
        self.tmin = tmin
        self.tmax = tmax

    # TODO:  Improve unit conversion
    # Now that is handled in thermodynamics
    def addProfile(
            self,
            pressure,
            temperature,
            dewPointTemperature,
            hPa=True,
            celsius=True,
            method=0,
            initialLevel=0,
            parcel=True,
            label=None,
            diagnostics=False,
            markers=True,
            useVirtual=True,
            **kwargs
    ):
        """
        Add a profile to the Skew-T axis

        .. _`matplotlib.legend`: http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.legend

        .. _ndarray: https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html

        .. _MaskedArray: https://docs.scipy.org/doc/numpy/reference/maskedarray.html

        The input data can  be ndarray_ or MaskedArray_.
        If MaskedArray_ are used, the pressure levels with masked or invalid entries
        in the pressure and temperature arrays are removed.
        Masked values in the dewPointTemperature are not plotted.

        Parameters
        ----------
        pressure : ndarray_ or MaskedArray_
            Pressure levels of the sounding
        temperature : ndarray_ or MaskedArray_
            Temperature profile at pressure levels
        dewPointTemperature : : ndarray_ or MaskedArray_
            Dew point temperature at pressure levels
        hPa: bool, optional
            If is True, the pressure levels are in hPa. Otherwise Pascals is
            assumed
        celsius : bool, optional
            If is True, the temperatures units correspond to celsius degrees.
            Otherwise Kelvin degrees units are assumed
        parcel : bool, optional
            If True, the parcel analysis is carried out.
        method : int, optional
            Parcel analysis method used. Supported:

            * Most Unstable  : method=0
            * Single Parcel: method=1
            * Mixed Layer  : method=2 (Not implemented)
        initialLevel : int, optional
            Initial level (index) used to compute the parcel analysis.
            Levels below this value are ignored.

            For the Single Parcel analysis, this level correspond to the parcel used.
            By default, the initial level is 0 (surface).

            For the Most Unstable method, this value is ignored.
        label : str, optional
            Label assigned to the profile lines. If None, no label is used.
            When a line is labeled, the legend is added automatically to the plot.
            use keyword **loc** to control the position.
        diagnostics : bool, optional
            If True a text box is added to the upper right corner with some diagnostics
            from the parcel analysis.
        markers: bool, optional
            If True, the LCL, LFC and EL are marked in the plot with markers.
        useVirtual : bool, optional
            If True, in the parcel analysis, CAPE and CIN are computed used
            the virtual temperature. The temperatures plotted in the SkewT diagram
            will correspond to virtual temperature instead of temperature.
            If False, virtual temperatures corrections are neglected and the
            original temperature is plotted.
        twColor: str
            Color for the dew point temperature line.
        tdColor: str
            Color for the dry temperature line.
        tpColor: str
            Color for the parcel temperature line.
        Other Parameters
        ----------------
        loc : str
            Legend location. See `matplotlib.legend`_ for more details.
        kwargs : extra
            The remaining extra keyword arguments are passed to the plot function
            when plotting the temperature profiles.
        """

        mask = numpy.zeros_like(pressure, dtype=bool)

        if isinstance(pressure, masked_array):
            pressure[pressure.mask] = numpy.nan
            mask = logical_or(mask, getmaskarray(pressure))
        if isinstance(temperature, masked_array):
            temperature[temperature.mask] = numpy.nan
            mask = logical_or(mask, getmaskarray(temperature))
        if isinstance(dewPointTemperature, masked_array):
            dewPointMask = getmaskarray(dewPointTemperature)
            dewPointTemperature.data[dewPointMask] = numpy.nan
            dewPointTemperature = numpy.asarray(
                dewPointTemperature.data, dtype=numpy.float32
            )

        pressure = masked_invalid(pressure[~mask])
        temperature = masked_invalid(temperature[~mask])
        dewPointTemperature = masked_invalid(dewPointTemperature[~mask])

        # Convert temperatures and pressure to hPa
        if not hPa:
            pressure /= 100

        if not celsius:
            temperature -= degCtoK
            dewPointTemperature -= degCtoK

        if parcel:

            parcelAnalysisResult = parcelAnalysis(
                pressure,
                temperature,
                dewPointTemperature,
                hPa=True,
                celsius=True,
                useVirtual=1,
                initialLevel=initialLevel,
                method=method,
            )

            initialLevel = parcelAnalysisResult["initialLevel"]
            initialTemperature = temperature[initialLevel]

            pressureAtLCL = parcelAnalysisResult["pressureAtLCL"]
            temperatureAtLCL = parcelAnalysisResult["temperatureAtLCL"]

            pressureAtLFC = parcelAnalysisResult["pressureAtLFC"]
            temperatureAtLFC = parcelAnalysisResult["temperatureAtLFC"]

            pressureAtEL = parcelAnalysisResult["pressureAtEL"]
            temperatureAtEL = parcelAnalysisResult["temperatureAtEL"]

            parcelTemperature = liftParcel(
                initialTemperature,
                pressure,
                pressureAtLCL,
                initialLevel=initialLevel,
                hPa=True,
                celsius=True,
            )

            # Add LCL
            belowLCL = numpy.where(pressure > pressureAtLCL, True, False)
            newParcelTemperature = numpy.concatenate(
                (
                    parcelTemperature[belowLCL],
                    [temperatureAtLCL],
                    parcelTemperature[~belowLCL],
                )
            )
            newPressure = numpy.concatenate(
                (pressure[belowLCL], [pressureAtLCL], pressure[~belowLCL])
            )

            # Add EL
            belowEL = numpy.where(newPressure > pressureAtEL, True, False)
            newParcelTemperature = numpy.concatenate(
                (
                    newParcelTemperature[belowEL],
                    [temperatureAtEL],
                    newParcelTemperature[~belowEL],
                )
            )
            newPressure = numpy.concatenate(
                (newPressure[belowEL], [pressureAtEL], newPressure[~belowEL])
            )

            belowLFC = numpy.where(newPressure > pressureAtLFC, True, False)
            newParcelTemperature = numpy.concatenate(
                (
                    newParcelTemperature[belowLFC],
                    [temperatureAtLFC],
                    newParcelTemperature[~belowLFC],
                )
            )
            newPressure = numpy.concatenate(
                (newPressure[belowLFC], [pressureAtLFC], newPressure[~belowLFC])
            )

            newTemperature = numpy.interp(
                newPressure, pressure[::-1], temperature[::-1]
            )
            newParcelTemperature = masked_invalid(newParcelTemperature)
            if useVirtual:
                newDewPointTemperature = numpy.interp(
                    newPressure, pressure[::-1], dewPointTemperature[::-1]
                )
                newTemperature = virtualTemp3(
                    newTemperature, newDewPointTemperature, newPressure * 100
                )

                belowLCL = newPressure >= pressureAtLCL

                newParcelTemperature[belowLCL] = virtualTemp3(
                    newParcelTemperature[belowLCL],
                    dewPointTemperature[initialLevel],
                    newPressure[belowLCL] * 100,
                )
                aboveLCL = newPressure < pressureAtLCL

                newParcelTemperature[aboveLCL] = virtualTemp4(
                    newParcelTemperature[aboveLCL], newPressure[aboveLCL] * 100
                )

        else:
            newTemperature = temperature
            newPressure = pressure

        kwargs["zorder"] = kwargs.pop("zorder", 5)
        kwargs["linewidth"] = kwargs.pop("linewidth", 2.0)
        kwargs["linestyle"] = kwargs.pop("linestyle", "-")
        loc = kwargs.pop("loc", "best")
        twColor = kwargs.pop("twColor", "b")
        tdColor = kwargs.pop("tdColor", "r")
        tpColor = kwargs.pop("tpColor", "k")

        (temperatureLine,) = self.plot(newTemperature, newPressure, tdColor, **kwargs)

        self.plot(dewPointTemperature, pressure, twColor, **kwargs)

        if label is not None:
            self.linesHandlers.append(temperatureLine)
            self.linesLabels.append(label)
            self.legend(self.linesHandlers, self.linesLabels, loc=loc)

        if parcel:
            self.plot(
                newParcelTemperature, newPressure, tpColor, **kwargs
            )

            if parcelAnalysisResult["CAPE"] > 0:
                # Positive Buoyancy
                cond1 = (newPressure <= pressureAtLFC) * (newPressure >= pressureAtEL)

                self.fill_betweenx(
                    newPressure,
                    newParcelTemperature,
                    newTemperature,
                    where=cond1,
                    color="#ff0009",
                    alpha=0.4,
                    zorder=10,
                )

                # Negative Buoyancy

                validMask = ~getmaskarray(masked_invalid(newParcelTemperature))

                cond2 = (newParcelTemperature[validMask] <= newTemperature[validMask]
                         ) * (newPressure[validMask] >= pressureAtLFC)

                self.fill_betweenx(
                    newPressure[validMask],
                    newParcelTemperature[validMask],
                    newTemperature[validMask],
                    where=cond2,
                    color="#045cff",
                    alpha=0.4,
                    zorder=10,
                )

            if markers:
                if useVirtual:
                    temperatureAtLCL = virtualTemp4(
                        temperatureAtLCL, pressureAtLCL * 100
                    )
                    temperatureAtLFC = virtualTemp4(
                        temperatureAtLFC, pressureAtLFC * 100
                    )
                    temperatureAtEL = virtualTemp4(temperatureAtEL, pressureAtEL * 100)

                self.plot(temperatureAtLCL, pressureAtLCL, ls="", marker="o", color="r")
                self.plot(temperatureAtLFC, pressureAtLFC, ls="", marker="o", color="g")
                self.plot(temperatureAtEL, pressureAtEL, ls="", marker="o", color="k")

            if diagnostics:
                # Add text to sounding
                dtext = "Diagnostics:\n"
                # dtext ="Parcel: %s\n"%ptype.upper()

                # dtext+="P  :%6.1fhPa\n"%startp
                # dtext+="T :  %4.1fC\n"%startt
                # dtext+="Td :  %4.1fC\n"%startdp
                dtext += "-------------\n"
                dtext += "P_par:%6.1fhPa\n" % (pressure[initialLevel])
                dtext += "P_LCL:%6.1fhPa\n" % (parcelAnalysisResult["pressureAtLCL"])
                dtext += "T_LCL:%4.1fC\n" % (parcelAnalysisResult["temperatureAtLCL"])
                dtext += "P_LFC:%6.1fhPa\n" % (parcelAnalysisResult["pressureAtLFC"])
                dtext += "T_LFC:%4.1fC\n" % (parcelAnalysisResult["temperatureAtLFC"])
                dtext += "P_EL:%6.1fhPa\n" % (parcelAnalysisResult["pressureAtEL"])
                dtext += "T_EL:%4.1fC\n" % (parcelAnalysisResult["temperatureAtEL"])
                dtext += "CAPE:%.1f\n" % parcelAnalysisResult["CAPE"]
                dtext += "CIN: %.1fJ" % parcelAnalysisResult["CIN"]

                axesBox = self.get_position().get_points()

                self.figure.text(
                    axesBox[1, 0],
                    axesBox[1, 1],
                    dtext,
                    fontname="monospace",
                    backgroundcolor="white",
                    zorder=10,
                    verticalalignment="top",
                    horizontalalignment="right",
                    multialignment="left",
                )

    def other_housekeeping(self, mixratio=array([])):
        """
        Set the Skew-T diagram properties:

        * y axis grid
        * x and y axis ticks
        * x axis label
        * x and y axis limits (see :py:meth:setLimits)
        """

        self.yaxis.grid(True, ls="-", color="y", lw=0.5)
        majorLocatorDegC = MultipleLocator(10)

        self.xaxis.grid(True, color="y", lw=0.5, ls="-")

        # self.set_ylabel('Pressure (hPa)')
        self.set_xlabel("Temperature (C)")
        yticks = linspace(100, 1000, 10)
        if self.pmin < 100:
            yticks = concatenate((array([50, 20, 10]), yticks))

        self.set_yticks(yticks, minor=False)
        self.set_yticks([], minor=True)

        yticks = linspace(100, 1000, 10)

        self.yaxis.set_major_formatter(ScalarFormatter())

        self.set_xlim(self.tmin, self.tmax)
        self.set_ylim(self.pmax, self.pmin)

        self.xaxis.set_major_locator(majorLocatorDegC)
        # self.spines['right'].set_visible(False)
        # self.get_yaxis().set_tick_params(which="minor",left='off')
        # self.get_xaxis().set_tick_params(which="both", size=0)

    def add_dry_adiabats(self, T0, P, do_labels=True, **kwargs):
        """ Add the dry adiabats to the Axis"""
        # Taken from Thomas Chubb's add_dry_adiabats method

        P0 = 1000.0
        T = array([(st + degCtoK) * (P / P0) ** (Rs_da / Cp_da) - degCtoK for st in T0])
        labelt = [(st + degCtoK) * 1 ** (Rs_da / Cp_da) for st in T0]

        # gets a pressure level about 2/4 the way up the plot...

        # pp = 10 ** (log10(self.pmin ** .2 * self.pmax ** .8))
        pp = 400  # at 400hPa

        xi = where(abs(P - pp) - abs(P - pp).min() < 1e-6)[0][0]

        ndec = log10(self.pmax / pp) / log10(self.pmax / self.pmin)
        tran = self.tmax - self.tmin
        tminl = self.tmin - tran * ndec
        tmaxl = self.tmax - tran * ndec

        if "color" in kwargs:
            col = kwargs["color"]
        else:
            col = "k"
        for tt, ll in zip(T, labelt):
            self.plot(tt, P, **kwargs)
            if do_labels:
                if tt[xi] > tmaxl - 2:
                    continue
                if tt[xi] < tminl + 2:
                    continue
                self.text(
                    tt[xi],
                    P[xi] + 10,
                    "%d" % (ll),
                    fontsize=8,
                    ha="center",
                    va="bottom",
                    rotation=-30,
                    color=col,
                    bbox={"facecolor": "w", "edgecolor": "w"},
                )
        return T

    def add_moist_adiabats(self, T0, P0, do_labels=True, nsteps=500, **kwargs):
        """
        Add the Moist Adiabats to the Axis
        """
        # Taken from Thomas Chubb's add_moist_adiabats method and modified to use the
        # cython moistAscent function

        T = zeros((len(T0), nsteps))

        for index, temp in enumerate(T0):
            P, T[index, :] = moistAscent(P0, temp, celsius=True, hPa=True)

        # gets a pressure level about 3/4 the way up the plot...
        pp = 10 ** (log10(self.pmin ** 0.75 * self.pmax ** 0.25))
        xi = where(abs(P - pp) - abs(P - pp).min() < 1e-6)[0][0]

        ndec = log10(self.pmax / pp) / log10(self.pmax / self.pmin)
        tran = self.tmax - self.tmin
        tminl = self.tmin - tran * ndec
        tmaxl = self.tmax - tran * ndec

        if "color" in kwargs:
            col = kwargs["color"]
        else:
            col = "k"

        for index in range(len(T0)):
            tt = T[index, :]
            self.plot(tt, P, **kwargs)

            # if (tt[-1]>-60) and (tt[-1]<-10):
            if do_labels:
                if tt[xi] > tmaxl - 2:
                    continue
                if tt[xi] < tminl + 2:
                    continue
                self.text(
                    tt[xi],
                    P[xi],
                    "%d" % (tt[0] + degCtoK),
                    ha="center",
                    va="bottom",
                    fontsize=8,
                    bbox={"facecolor": "w", "edgecolor": "w"},
                    color=col,
                )

    def add_mixratio_isopleths(self, w, P, do_labels=True, **kwargs):
        """ Add the vapor mixing ration isopleths to the Axis"""
        # Taken from Thomas Chubb's add_dry_adiabats method

        if len(P) > 0:
            e = array([P * ww / (0.622 + ww) for ww in w])
            T = 243.5 / (17.67 / log(e / 6.112) - 1)
            if "color" in kwargs:
                col = kwargs["color"]
            else:
                col = "k"

            pp = 700.0
            xi = where(abs(P - pp) - abs(P - pp).min() < 1e-6)[0][0]

            ndec = log10(self.pmax / pp) / log10(self.pmax / self.pmin)
            tran = self.tmax - self.tmin
            tminl = self.tmin - tran * ndec
            tmaxl = self.tmax - tran * ndec

            for tt, mr in zip(T, w):
                self.plot(tt, P.flatten(), **kwargs)
                if do_labels:
                    if tt[xi] > tmaxl - 2:
                        continue
                    if tt[xi] < tminl + 2:
                        continue
                    if mr * 1000 < 0.1:
                        fmt = "%4.2f"
                    elif mr * 1000 <= 1.0:
                        fmt = "%4.1f"
                    else:
                        fmt = "%d"
                    self.text(
                        tt[-1],
                        P[-1],
                        fmt % (mr * 1000),
                        color=col,
                        fontsize=8,
                        ha="center",
                        va="bottom",
                        bbox={"facecolor": "w", "edgecolor": "w"},
                    )


# Now register the projection with matplotlib so the user can select it.
register_projection(SkewXAxes)


def figure(*args, **kwargs):
    """
    Creates a new figure along with the Figure Managers.

    **show** method is now part of the Figure class.

    Important!: When an instance is created, the Figure doesn't
    contains any axes.
    They must be created a posteriori using add_subplot_ for example.

    .. _add_subplot: https://matplotlib.org/api/figure_api.html#matplotlib.figure.Figure.add_subplot

    Parameters
    ----------

    num : integer or string, optional, default: none
        If not provided, a new figure will be created, and the figure number
        will be incremented. The figure objects holds this number in a `number`
        attribute.
        If num is provided, and a figure with this id already exists, make
        it active, and returns a reference to it. If this figure does not
        exists, create it and returns it.
        If num is a string, the window title will be set to this figure's
        `num`.

    figsize : tuple of integers, optional, default: None
        width, height in inches. If not provided, defaults to rc
        figure.figsize.

    dpi : integer, optional, default: None
        resolution of the figure. If not provided, defaults to rc figure.dpi.

    facecolor :
        the background color. If not provided, defaults to rc figure.facecolor

    edgecolor :
        the border color. If not provided, defaults to rc figure.edgecolor


    Returns
    -------
    figure : Figure
        The Figure instance returned will also be passed to new_figure_manager
        in the backends, which allows to hook custom Figure classes into the
        pylab interface. Additional kwargs will be passed to the figure init
        function.

    Notes
    -----
    If you are creating many figures, make sure you explicitly call "close"
    on the figures you are not using, because this will enable pylab
    to properly clean up the memory.

    rcParams defines the default values, which can be modified in the
    matplotlibrc file

    """
    my_figure = pyplot.figure(*args, **kwargs)

    def show_plot(*args, **kwargs):
        my_figure.canvas.draw()
        my_figure.canvas.flush_events()
        pyplot.show(*args, **kwargs)

    my_figure.show_plot = show_plot

    def save_fig(*args, **kwargs):
        my_figure.canvas.draw()
        my_figure.canvas.flush_events()
        pyplot.savefig(*args, **kwargs)

    my_figure.save_fig = save_fig

    return my_figure
