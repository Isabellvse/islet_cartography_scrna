from pybiomart import Server

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
