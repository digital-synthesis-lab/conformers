import click
from conformers.cli.density import density
from conformers.cli.process_list import process_list


class ConfGroup(click.Group):
    pass


@click.command(cls=ConfGroup)
def conf():
    """Command line interface for conformers"""


conf.add_command(density)
conf.add_command(process_list)

if __name__ == "__main__":
    conf()
