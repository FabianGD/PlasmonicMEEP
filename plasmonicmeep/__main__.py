"""Main entrypoint for PlasmonicMEEP restructuring"""

import argparse
import sys


def _do_fdtd(args):
    print("FDTD")
    ...


def _do_enhancement(args):
    print("ENH")
    ...


def _do_visualize(args):
    print("VIS")
    ...


def argparsing(args):
    parser = argparse.ArgumentParser(description="Set of scripts for calculation of plasmonic field enhancement, powered by MEEP.")

    parser.add_argument("-v", "--verbose", action="count")
    subparsers = parser.add_subparsers()

    fdtd_parser = subparsers.add_parser("fdtd")
    fdtd_parser.set_defaults(func=_do_fdtd)

    enh_parser = subparsers.add_parser("enh")
    enh_parser.set_defaults(func=_do_enhancement)

    vis_parser = subparsers.add_parser("vis")
    vis_parser.set_defaults(func=_do_visualize)

    return parser.parse_args(args)


if __name__ == "__main__":

    args = sys.argv[1:]
    p_args = argparsing(args)
    p_args.func(p_args)