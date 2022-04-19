def calc_extent(self, plot_f=None, plot_t=None, MJD_time=False):
    """ Setup plotting edges.
    """

    plot_f_begin = plot_f[0]
    plot_f_end = plot_f[-1] + (plot_f[1] - plot_f[0])

    plot_t_begin = self.timestamps[0]
    plot_t_end = self.timestamps[-1] + (self.timestamps[1] - self.timestamps[0])

    if MJD_time:
        extent = (plot_f_begin, plot_f_end, plot_t_begin, plot_t_end)
    else:
        extent = (plot_f_begin, plot_f_end, 0.0, (plot_t_end - plot_t_begin) * 24. * 60. * 60)

    return extent