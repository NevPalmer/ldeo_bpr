"""Test for arg_parser module."""

import ldeo_bpr as bpr

args = bpr.parse_arguments()

for arg in args.__dict__.items():
    print(arg)

print(__file__)
