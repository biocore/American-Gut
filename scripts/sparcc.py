# python 2
# conda install biom-format -c bioconda
# conda install matplotlib  # pysurvey setup.py is incomplete
# pip install -e "hg+https://bitbucket.org/yonatanf/pysurvey#egg=pysurvey"
# note: the pysurvey build in pypi does not work. unsure why.

import click

# this little gem is because pysurvey depends on an older version of scipy
# and attempts to import nanmedian from scipy.stats. This is a dirty hack.
import scipy.stats as ss
import scipy
ss.nanmedian = scipy.nanmedian

import pysurvey as ps
import biom
import pandas as pd


def _biom_to_pysurvey_mat(table):
    """Convert a BIOM table to a compatible pysurvey DataFrame

    Parameters
    ----------
    table : biom.Table
        The BIOM table

    Returns
    -------
    DataFrame
        A pandas DataFrame representing the BIOM table where the rows are
        samples and the columns are observations. The sample identifiers are
        stripped as they are not used, and any BIOM metadata are ignored
    """
    table = biom.load_table(table)
    mat = table.matrix_data.toarray().T
    return pd.DataFrame(mat, columns=table.ids(axis='observation'),
                        index=table.ids())


def _format_edges(fp, sparcc_z, threshold):
    """Generate edges

    Parameters
    ----------
    fp : open file
        Where to write the results
    sparcc_z : pd.DataFrame
        A pandas DataFrame containing the correlations
    threshold : float
        A minimum absolute correlation threshold value
    """
    fp.write('Feature1\tFeature2\tRho\tRho_pos_neg\n')
    for i, (idx, row) in enumerate(sparcc_z.iterrows()):

        # i+1 for upper triangle
        for col, v in zip(sparcc_z.columns[i+1:], row[i+1:]):
            if abs(v) >= threshold:
                fp.write('%s\t%s\t%f\t%d\n' % (idx, col, v, v > 0))


def _format_nodes(fp, table):
    """Format the nodes file

    Parameters
    ---------
    fp : open file
        Where to write the results too
    table : file path
        A file path to the BIOM table
    """
    header = ['Feature1', 'Fsum']

    table = biom.load_table(table)
    _, _, tmp_md = next(table.iter(axis='observation'))
    header.extend(sorted(tmp_md.keys()))

    fp.write("\t".join(header))
    fp.write('\n')

    for values, id_, md in table.iter(axis='observation'):
        line = [str(id_), str(values.sum())]
        for key in sorted(md.keys()):
            md_value = md[key]
            if isinstance(md_value, (list, tuple, set)):
                line.append(" ".join([str(v) for v in md_value]))
            else:
                line.append(str(md_value))
        fp.write('\t'.join(line))
        fp.write('\n')


@click.group()
def sparcc():
    pass


@sparcc.command()
@click.option('--table', required=True,
              type=click.Path(exists=True, dir_okay=False),
              help='The input BIOM table')
@click.option('--output', required=True,
              type=click.Path(exists=False),
              help='Filepath for the resulting correlation matrix')
def correlations(table, output):
    """Obtain SparCC correlations"""
    mat = _biom_to_pysurvey_mat(table)
    sparcc_z, cov = ps.basis_corr(mat)
    sparcc_z.to_csv(output, sep='\t', index_label='#OTU ID')


@sparcc.command()
@click.option('--correlations', required=True,
              type=click.Path(exists=True, dir_okay=False),
              help='The input correlation matrix')
@click.option('--output', required=True,
              type=str,
              help='Filename base for the resulting edge and node files')
@click.option('--threshold', required=False, default=0.35, type=float,
              help='The correlation threshold (applied to the absolute value)')
@click.option('--table', required=True,
              type=click.Path(exists=True),
              help='BIOM table with OTU metadata')
def network(correlations, output, threshold, table):
    """Construct network from correlations"""
    sparcc_z = pd.read_csv(correlations, sep='\t', index_col=0)

    with open(output + '.edges', 'w') as fp:
        _format_edges(fp, sparcc_z, threshold)

    with open(output + '.nodes', 'w') as fp:
        _format_nodes(fp, table)


if __name__ == '__main__':
    sparcc()
