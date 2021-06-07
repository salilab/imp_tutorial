#!/usr/bin/env python

"""Simple script to plot data from IMP/PMI stat files"""

import optparse
import IMP.pmi.output


class StatFileProcessor(object):
    def __init__(self, fname):
        self.fname = fname

    def _get_field_mapping(self, fh):
        """Parse the first line of the stat file to get all fields"""
        header = fh.readline()
        id_to_field_name = eval(header, {}, {})
        # Get reverse mapping
        field_name_to_id = {}
        for key, val in id_to_field_name.items():
            field_name_to_id[val] = key
        return id_to_field_name, field_name_to_id

    def get_names_values(self, fields):
        """Yield values for the given fields, for all frames in the file"""
        with open(self.fname) as fh:
            id_to_field_name, field_name_to_id = self._get_field_mapping(fh)
            fields = [self._map_field(x, id_to_field_name, field_name_to_id)
                      for x in fields]
            yield [f[0] for f in fields]
            ids = [f[1] for f in fields]
            for n, line in enumerate(fh.readlines()):
                data = eval(line, {}, {})
                data['Frame'] = n + 1
                yield [data[i] for i in ids]

    def _map_field(self, field_name, id_to_field_name, field_name_to_id):
        """Given an input field name, return name and id"""
        if field_name == 'Frame':
            return 'Frame', 'Frame'
        try:
            field_id = int(field_name)
            return id_to_field_name[field_id], field_id
        except ValueError:
            pass
        return field_name, field_name_to_id[field_name]


def parse_args():
    usage = """%prog [options] <stat file>

Simple script to plot data from IMP/PMI stat files.
"""
    parser = optparse.OptionParser(usage)
    parser.add_option("-x", metavar='FIELD', default="Frame",
                      help="Field to plot on the x axis. This should be a "
                           "column name, such as 'GaussianEMRestraint_None', "
                           "a column number, or 'Frame' to plot the frame "
                           "number.")
    parser.add_option("-y", metavar='FIELD', default="Frame",
                      help="Field to plot on the y axis. See help for -x.")
    parser.add_option("-b", metavar='FRAME', default=1, dest='begin', type=int,
                      help="Begin from this frame (in case you want to skip "
                           "an equilibration phase")
    opts, args = parser.parse_args()
    if len(args) != 1:
        parser.error("Incorrect number of arguments")
    return opts, args[0]


def main():
    opts, statfile = parse_args()
    f = StatFileProcessor(statfile)
    outs = list(f.get_names_values([opts.x, opts.y]))
    xlabel, ylabel = outs[0]
    IMP.pmi.output.plot_xy_data([o[0] for o in outs[opts.begin:]],
                                [o[1] for o in outs[opts.begin:]],
                                out_fn="%s-%s" % (xlabel, ylabel),
                                xlabel=xlabel, ylabel=ylabel)


if __name__ == '__main__':
    main()
