#! /usr/bin/env python3

import click
import os
from CheckLengths import checkLengths


def globalCheckLengths(input:str):

    files = os.listdir(input)
    for file in files:
        local_input = os.path.join(input, file)
        checkLengths(local_input)


@click.command()
@click.argument('input', type=click.Path(exists=True))


def execute(input):

    globalCheckLengths(input)


if __name__ == '__main__':
    execute()
