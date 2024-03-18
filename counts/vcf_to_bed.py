# import vcfpy
from pysam import VariantFile
from io import StringIO
import os
import shutil
import click

@click.command()
@click.argument("vcf")
@click.option("-p", "--padding", default=1000, help="pad +/- from start ends")
@click.option("-o", "--out", default="beds", help="out dir")
def vcf_to_bed(vcf, padding, out):
    if os.path.isfile(vcf):
        files = [vcf]
    if os.path.isdir(vcf):
        files = [os.path.join(vcf, f) for f in os.listdir(vcf) if f.endswith(".vcf")]
    total_files = len(files)
    print(f"{total_files}: {files}")
    for num, file in enumerate(files):
        print(f"{num}/{total_files} ({round((num/total_files)*100, 2)}%) {file}", end="\r")
        f = VariantFile(file)
        out_file = StringIO()
        for v in f:
            if v.stop < v.pos:
                print(v)
                continue
            if v.stop - v.pos <= 1000:
                out_file.write(f"{v.chrom}\t{0 if v.pos-padding < 0 else v.pos-padding}\t{v.stop+padding}\n")
            else:
                out_file.write(f"{v.chrom}\t{0 if v.pos-padding < 0 else v.pos-padding}\t{v.pos+padding}\n")
                out_file.write(f"{v.chrom}\t{0 if v.stop-padding < 0 else v.stop-padding}\t{v.stop+padding}\n")
            if "CHR2_pos" in v.info:
                out_file.write(f"{v.info['CHR2']}\t{v.info['CHR2_POS']-padding}\t{v.info['CHR2_POS']+padding}\n")
        with open(f"{os.path.join(out, os.path.splitext(os.path.basename(file))[0])}.bed", mode="w") as of:
            out_file.seek(0)
            shutil.copyfileobj(out_file, of)
            of.close()

if __name__ == "__main__":
    vcf_to_bed()

