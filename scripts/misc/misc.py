import re
import matplotlib            # Base matplotlib functionality
import seaborn as sns
import matplotlib.pyplot as plt
from pybiomart import Server

def to_snake_case(name):
    # Convert to lowercase, replace spaces and parentheses with underscores, remove multiple underscores
    name = name.lower()
    name = re.sub(r"[^\w]+", "_", name)
    name = re.sub(r"_+", "_", name)
    return name.strip("_")

def set_my_theme():
    sns.set_style("white")  # similar to theme_classic background

    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42

    plt.rcParams.update({
        # text + labels
        "axes.titlesize": 8,
        "axes.titleweight": "normal",
        "axes.labelsize": 7,
        "xtick.labelsize": 6,
        "ytick.labelsize": 6,
        "legend.fontsize": 4,
        "legend.title_fontsize": 4,

        # line & tick colors
        "axes.edgecolor": "black",
        "xtick.color": "black",
        "ytick.color": "black",
        "text.color": "black",

        # spine & tick thickness
        "axes.linewidth": 0.7,
        "xtick.major.width": 0.7,
        "ytick.major.width": 0.7,

        # IMPORTANT: make ticks long & visible
        "xtick.major.size": 4,
        "ytick.major.size": 4,
        "xtick.bottom": True,
        "ytick.left": True,

        # background
        "axes.facecolor": "white",
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
    })


def mouse_2_human(df, mouse_col='mouse_symbol'):
    """
    Parameters:
    df (pd.DataFrame): DataFrame containing mouse genes
    mouse_col (str): Column name for mouse gene symbols

    Returns:
        pd.DataFrame: Original DataFrame with an additional column for human orthologs
    """
    # Connect to Ensembl
    server = Server(host='www.ensembl.org')
    mouse_dataset = server.marts['ENSEMBL_MART_ENSEMBL'].datasets['mmusculus_gene_ensembl']
    
    # Get unique mouse genes
    mouse_genes = df[mouse_col].unique().tolist()
    
    # Query Ensembl for human orthologs
    gene_map = mouse_dataset.query(
        attributes=[
            'external_gene_name',                       # mouse gene symbol
            'hsapiens_homolog_associated_gene_name',   # human ortholog symbol
            'hsapiens_homolog_orthology_type',         # 1:1, 1:many, many:many
            'hsapiens_homolog_perc_id'                 # protein sequence similarity %
        ]
    )

    # keep only genes we are interested in
    gene_map.query("`Gene name` in @mouse_genes")
    
    # Remove missing human orthologs
    gene_map = gene_map.query("`Human gene name`.notna() and `Human gene name` != '' and `Human gene name` != 'nan'")
    
    # Prefer 1:1 orthologs, then highest percent identity
    gene_map = (gene_map
                .sort_values(by=['Gene name',
                                 'Human homology type',
                                 '%id. target Human gene identical to query gene'], 
                             ascending=[True, False, False])
                .groupby('Gene name', as_index=False)
                .first())
    
    # Merge back with original df
    df_mapped = df.merge(
        gene_map[['Gene name', 'Human gene name']]
        .rename(columns={'Human gene name': 'human_symbol'}),
        left_on=mouse_col,
        right_on='Gene name',
        how='left'
    ).drop(columns=['Gene name'])
    
    # Remove duplicates in human genes
    df_mapped = df_mapped.drop_duplicates(subset=['human_symbol'])

    # Remove missing human symbols
    df_mapped = df_mapped.query("human_symbol.notna() and human_symbol != '' and human_symbol != 'nan'")
    
    return df_mapped
