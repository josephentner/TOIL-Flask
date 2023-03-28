from flask import Flask, request
import json
import pandas as pd
import xenaPython as xena

cohort = "TCGA TARGET GTEx"
host = xena.PUBLIC_HUBS["toilHub"]

def get_field(dataset, samples, probe):
    """
    Returns all values for a specific field for a list of samples

    Parameters:
        dataset: name of Xena dataset to query (string)
        samples: list of sample ID"s (list)
        probe: name of field to get values for (string)
    
    Returns:
        field_values: list of values for specified samples (list)
    """

    # get dictionary to map result code to actual results
    field_codes = xena.field_codes(host, dataset, [probe])[0]["code"]
    codes_dict = dict(enumerate(field_codes.split("\t")))
    
    # get field value codes and map to results
    raw_values = xena.dataset_fetch(host, dataset, samples, [probe])[0]
    field_values = list(map(codes_dict.get, raw_values))

    return field_values

def get_data(gene, gene2, diseases):
    """
    Gets expression results and relevant metadata for all samples for a specific gene

    Parameters:
        gene: name of gene to query (string)
    
    Returns: 
        data: dataframe that contains all samples, expression results, and metadata (pd.DataFrame)
    """

    # get expression results
    dataset = "TcgaTargetGtex_RSEM_Hugo_norm_count"
    samples = xena.dataset_samples(host, dataset, None)
    expression = xena.dataset_gene_probes_values(host, dataset, samples, gene)[1][0]
    
    # get relevant metadata (sample type and disease/tissue)
    dataset = "TcgaTargetGTEX_phenotype.txt"
    disease_tissue = get_field(dataset, samples, "primary disease or tissue")

    # create DataFrame, clean data results, and store data
    data = pd.DataFrame()
    data["Sample"] = samples
    data["Study"] = data["Sample"].apply(lambda x: x.split("-")[0])
    data["Expression"] = expression
    data["Disease/Tissue"] = disease_tissue
    data["Disease/Tissue"] = data["Disease/Tissue"].apply(lambda x : x.split(" - ")[0] if x is not None else x )
    
    if gene2 != '':
        expression2 = xena.data_gene_probes_values(host, dataset, samples, gene2)[1][0]
        data["Expression2"] = expression2

    # filter data 
    graph_data = data[(data["Disease/Tissue"].isin(list(diseases)))]
    
    # normal_data = data[(data["Study"] == "GTEX") & (data["Sample Type"].isin(sample_types)) & data["Disease/Tissue"].isin(tissues)]
    # cancer_data = data[(data["Disease/Tissue"].isin(list(diseases))) & data["Sample Type"].isin(sample_types)]
    # graph_data = pd.concat([cancer_data, normal_data], axis=0)
    
    return graph_data


def get_genes():
    """
    Grabs all genes from the dataset

    Parameters: None

    Returns: 
        genes: list of all genes (list)
    """
    
    dataset = "TcgaTargetGtex_RSEM_Hugo_norm_count"
    genes = xena.dataset_field(host, dataset)

    return genes

def get_diseases():
    """
    Grabs all diseases from the dataset

    Parameters: None

    Returns: TCGA diseases and TARGET diseases (list)
    """
    
    host = xena.PUBLIC_HUBS["toilHub"]
    dataset = "TcgaTargetGTEX_phenotype.txt"
    probe = "primary disease or tissue"

    # create dictionary that stores all diseases as values
    field_codes = xena.field_codes(host, dataset, [probe])[0]["code"]
    codes_dict = dict(enumerate(field_codes.split("\t")))

    # keep only the diseases that are in the TCGA and TARGET datasets
    tcga_diseases = list(codes_dict.values())[:33]
    target_diseases = list(codes_dict.values())[88:]

    return tcga_diseases + target_diseases


app = Flask(__name__)

@app.route('/')
def layout():
    return 'Hello World'

@app.route('/data', methods=['GET'])
def send_data():
    diseases = request.args.get('disease', default=['Breast Invasive Carcinoma'], type=str)
    gene = request.args.get('gene', default='ERBB2', type=str)
    gene2 = request.args.get('gene2', default='', type=str)
    data = get_data(gene, gene2, diseases)
    return data.to_json()

@app.route('/genes', methods=['GET'])
def send_genes():
    genes = get_genes()
    return json.dumps(genes)

@app.route('/diseases', methods=['GET'])
def send_diseases():
    diseases = get_diseases()
    return json.dumps(diseases)

if __name__ == '__main__':
    app.run(debug = True, host = '0.0.0.0', port = 8000)

