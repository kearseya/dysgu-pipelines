import click

import fetch_and_assemble
import map_reads_to_contigs
import analyse_contigs
import remove_common_sites
import collect_mapping_info
import plot_contigs

@click.group(name="pipe", invoke_without_command=True)
@click.option("-b", "--bam-glob", default="/mnt/breast/*.*am")
@click.option("-c", "--chain-file", default="../chain_link/chain_out_1/t20_p0.05-0.35/found_chains_p0.2_t20/all_svs.unique.chains.csv")
@click.option("-r", "--ref", defualt="/mnt/scwc0010/hg38.fa")
@click.option("-p", "--pad", default=500)
@click.pass_context
def pipe(ctx, bam_glob, chain_filea, ref, pad):
    ctx.ensure_object(dict)
    ctx.obj["bam_glob"] = bam_glob
    ctx.obj["chain_file"] = chain_file
    ctx.obj["ref"] = ref
    ctx.obj["pad"] = pad
    if ctx.invoked_subcommand is None:
        if bam_glob == None:
            click.echo(ctx.get_help())
            ctx.exit()
            pass
        else:
            click.echo("Running pipeline")
            ctx.invoke(fetch_and_assemble_pipe)
            ctx.invoke(map_reads_to_contigs_pipe)
            ctx.invoke(analyse_contigs_pipe)
            ctx.invoke(remove_common_sites_pipe)
    else:
        pass

@pipe.command(name="fetch")
@click.pass_context
def fetch_and_assemble_pipe(ctx):
    """ 1st fetch and assemble contigs """
    print("Fetching and assembling reads...")
    fetch_and_assemble.main(ctx.obj["bam_glob"], ctx.obj["chain_file"], ctx.obj["ref"], ctx.obj["pad"])
    pass

@pipe.command(name="map")
@click.pass_context
def map_reads_to_contigs_pipe(ctx):
    """ 2nd map reads to contigs """
    print("Mapping reads to contigs...")
    map_reads_to_contigs.main("regions_all")
    pass

@pipe.command(name="analyse")
@click.pass_context
def analyse_contigs_pipe(ctx):
    """ 3rd analyse contigs """
    print("Analysing contigs...")
    analyse_contigs.main("regions_all", ref=ctx.obj["ref"])
    pass

@pipe.command(name="remove")
@click.pass_context
def remove_common_sites_pipe(ctx):
    """ 4th remove common sites """
    print("Removing common sites...")
    remove_common_sites.load_data()
    pass

@pipe.command(name="collect")
@click.pass_context
def collect_mapping_info(ctx):
    """ 5th collet read map """
    print("Collecting read mapping...")
    collect_mapping_info.mapping_info("out_data/all.bwa_dodi.bam", "out_data/contigs.bed", ctx.obj["chain_file"], ctx.obj["pad"])

@pipe.command(name="plot")
@click.pass_context
def plot_all_contigs(ctx):
    """ 6th plot contigs """
    print("Plotting contigs...")
    plot_contigs.load_data("out_data/filtered_contigs.bed", ctx.obj["ref"])

if __name__ == "__main__":
    pipe()
