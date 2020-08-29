###### CODE SETUP 

## First get all the functions set up
import pandas as pd
import requests
import difflib
import math

## Load BTE
from biothings_explorer.user_query_dispatcher import FindConnection
from biothings_explorer.hint import Hint
ht = Hint()

# all intermediate node types
# node_type_list = (['Gene', 'SequenceVariant', 'ChemicalSubstance', 'Disease', 
#                 'MolecularActivity', 'BiologicalProcess', 'CellularComponent', 
#                 'Pathway', 'AnatomicalEntity', 'PhenotypicFeature'])

# node_type_list = (['PhenotypicFeature'])

# Check for list of intermediate node types node type using predict funciton
def predict_many(input_object, intermediate_node_list, output_type):
    df_list = []
    for inter in intermediate_node_list:
        try: 
            print("Intermediate Node type running:")
            print(inter)
            fc = FindConnection(input_obj=input_object, output_obj=output_type, intermediate_nodes=[inter])
            fc.connect(verbose=False)
            df = fc.display_table_view()
            rows = df.shape[0]
            if(rows > 0):
                df_list.append(df)
        except:
            print("FAILED")
    if(len(df_list) > 0):
        return pd.concat(df_list)
    else:
        return None



def get_disease_to_gene_results(disease_input):
    disease_to_gene_results = {}

    #directly related
    fc = FindConnection(input_obj=disease_input, output_obj='Gene', intermediate_nodes=None)
    fc.connect(verbose=False)
    disease_to_genes = fc.display_table_view()

    disease_to_genes = disease_to_genes[~disease_to_genes['output_id'].str.contains('UMLS')]

    # keep track of number of occurrences from direct disease -> gene connection
    print("running disease -> gene")
    i = list(disease_to_genes["output_name"])
    d = {x:i.count(x) for x in i}
    sorted_disease_to_genes = {k: v for k, v in sorted(d.items(), key=lambda item: item[1])}
    disease_to_gene_results["sorted_disease_to_genes"] = sorted_disease_to_genes
    # print("occurences of genes directly related to genes")
    # print(disease_to_gene_results["sorted_disease_to_genes"])

    one_step_genes_pub_counts = {}
    for index, row in disease_to_genes.iterrows():
        current_pubcount = 0
        if(row["pred1_pubmed"] != None):
            current_pubcount = current_pubcount + row["pred1_pubmed"].count(",") + 1
        if row["output_name"] in one_step_genes_pub_counts:
            one_step_genes_pub_counts[row["output_name"]] = one_step_genes_pub_counts[row["output_name"]] + current_pubcount
        else: 
            one_step_genes_pub_counts[row["output_name"]] = current_pubcount
    disease_to_gene_results["one_step_genes_pub_counts"] =  one_step_genes_pub_counts

    disease_to_genes_list = list(reversed(list(sorted_disease_to_genes.keys())))
    disease_to_gene_results["disease_to_genes_list"] = disease_to_genes_list

    return(disease_to_gene_results)



def get_disease_to_node_to_gene_results(disease_all_nodes_genes,max_two_step_gene_count,symptom_list,symptoms_hpids):
    disease_to_node_to_gene_results = {}

    print("finding intermediate nodes that are symptoms")
    indices_with_symptom_as_intermediate = []
    go_dict = {}
    for index, row in disease_all_nodes_genes.iterrows():
        if row["node1_type"] == 'Disease':
            if row["node1_name"].lower() in symptom_list:
                indices_with_symptom_as_intermediate.append(index)
        elif row["node1_type"] == 'BiologicalProcess':
            if row["node1_name"].lower() in symptom_list:
                indices_with_symptom_as_intermediate.append(index)
            elif("go:" in row["node1_name"]):
                go_ID = row["node1_name"].split(':')[1]
                if go_ID in go_dict: 
                    if go_dict[go_ID] == True:
                        indices_with_symptom_as_intermediate.append(index)
                else:
                    r = requests.get('https://biothings.ncats.io/go_bp/geneset/GO%3A' + go_ID)
                    res = r.json()
                    if('name' in res):
                        if res['name'] in symptom_list: 
                            indices_with_symptom_as_intermediate.append(index)
                            go_dict[go_ID] = True
                        else:
                            go_dict[go_ID] = False
        elif row["node1_type"] == 'PhenotypicFeature':
            if row["node1_name"].lower() in symptom_list:
                indices_with_symptom_as_intermediate.append(index)
            elif("HP" in row["node1_name"]):
                if row["node1_name"] in symptoms_hpids:
                    indices_with_symptom_as_intermediate.append(index)

    # print("indices")
    # print(indices_with_symptom_as_intermediate)
    print("removing symptom intermediates")
    disease_all_nodes_genes = disease_all_nodes_genes.drop(disease_all_nodes_genes.index[indices_with_symptom_as_intermediate])

    print("getting gene counts from " + str(len(list(disease_all_nodes_genes["output_name"]))) + " gene entries" )
    i = list(disease_all_nodes_genes["output_name"])
    # d = {x:i.count(x) for x in i}
    d = {}
    for x in i: 
        if x in d: 
            d[x] = d[x] + 1
        else:
            d[x] = 1
    print("sorting counts dictionary")
    sorted_disease_to_all_nodes_to_genes = {k: v for k, v in sorted(d.items(), key=lambda item: item[1])}

    print("top genes occurrence counts: ")
    for x in list(reversed(list(sorted_disease_to_all_nodes_to_genes)))[0:max_two_step_gene_count]:
        print(str(x) + ": " + str(sorted_disease_to_all_nodes_to_genes[x]))
    
    top_related_genes_to_disease = list(reversed(list(sorted_disease_to_all_nodes_to_genes)))[0:max_two_step_gene_count]

    disease_to_node_to_gene_results["top_related_genes_to_disease"] = top_related_genes_to_disease
    disease_to_node_to_gene_results["sorted_disease_to_all_nodes_to_genes"] = sorted_disease_to_all_nodes_to_genes

    # keep track of pubication counts for genes in two-step disease -> intermediate node -> gene
    print("getting publicaiton counts")
    top_two_step_genes_pub_counts = {}
    for index, row in disease_all_nodes_genes.iterrows():
        if row["output_name"] in top_related_genes_to_disease:
            current_pubcount = 0
            if(row["pred1_pubmed"] != None):
                current_pubcount = current_pubcount + str(row["pred1_pubmed"]).count(",") + 1
            if(row["pred2_pubmed"] != None):
                current_pubcount = current_pubcount + str(row["pred2_pubmed"]).count(",") + 1
            if row["output_name"] in top_two_step_genes_pub_counts:
                top_two_step_genes_pub_counts[row["output_name"]] = top_two_step_genes_pub_counts[row["output_name"]] + current_pubcount
            else: 
                top_two_step_genes_pub_counts[row["output_name"]] = current_pubcount

    disease_to_node_to_gene_results["top_two_step_genes_pub_counts"] =  top_two_step_genes_pub_counts
    
    return(disease_to_node_to_gene_results)


def get_disease_symptoms(disease_name):
    request_url = 'http://mydisease.info/v1/query?q=hpo.disease_name:"' + disease_name + '"&fields=hpo'
    print("request url: " + request_url)
    r = requests.get(request_url)
    res = r.json()
    # print(res)
    result_number = 0
    disease_info = res['hits'][result_number]
    print("disease symptoms for:")
    print(disease_info['hpo']['disease_name'])
    symptoms = []
    hp_ids = []
    hp_symptom_dict = {}
    for x in disease_info['hpo']['phenotype_related_to_disease']:
        if(x["aspect"].lower() == 'p'):
            try:
                # if('frequency' in x):
                r1 = requests.get('https://biothings.ncats.io/hpo/phenotype/' + x['frequency'])
                res1 = r1.json()
                frequency = True
            except:
                frequency = False; 
                # print('no frequency')
            r = requests.get('https://biothings.ncats.io/hpo/phenotype/' + x['hpo_id'])
            res = r.json()
            if(('_id' in res) & ('name' in res)):
                symptoms.append(res['name'].lower())
                hp_ids.append(res['_id'])
                hp_symptom_dict[res['_id']] = {
                    'names' : [res['name'].lower()],
                    'frequency' : res1['name'] if frequency else 'Unknown',
                    'edges_out_count' : 0
                }
            if('synonym' in res):
                if('exact' in res['synonym']):
                    for name in res['synonym']['exact']:
                        if name not in symptoms: 
                            symptoms.append(name)
                            hp_symptom_dict[res['_id']]['names'].append(name)

    print(symptoms)
    return([symptoms,hp_ids,hp_symptom_dict])


def get_symtpom_prevalence(hp_symptom_dict, disease_name):
    for key in hp_symptom_dict:
        print(key)
        edges_out_count = 0
        # print("name: " + str(hp_symptom_dict[key]))
        UMLS = ''
        for y in ['PhenotypicFeature','Disease','BiologicalProcess']:
            if y == 'PhenotypicFeature':
                a = ht.query(key)[y]
                if len(a) > 0: 
                    b = a[0]
                    if 'UMLS' in b: 
                        UMLS = b['UMLS']
                    try: 
                        fc = FindConnection(input_obj=b, output_obj='Gene', intermediate_nodes=None)
                        fc.connect(verbose=False)
                        df = fc.display_table_view()
                        print('gene')
                        print(df.shape)
                        # print("phen")
                        # print(hp_symptom_dict[key])
                        # print(df.shape[0])
                        if(df.shape[0] > 0):
                            print("OKKKk")
                            df = df[df["output_name"] != disease_name]
                            edges_out_count = edges_out_count + df.shape[0]
                            print(edges_out_count)
                    except: 
                        print("Nope")
                    try: 
                        fc = FindConnection(input_obj=b, output_obj='Disease', intermediate_nodes=None)
                        fc.connect(verbose=False)
                        df = fc.display_table_view()
                        print(df.shape)
                        # print("phen")
                        # print(hp_symptom_dict[key])
                        # print(df.shape[0])
                        if(df.shape[0] > 0):
                            print("ok edge phen to dis")
                            df = df[df["output_name"] != disease_name]
                            edges_out_count = edges_out_count + df.shape[0]
                            print(edges_out_count)
                    except: 
                        print("Nope")
            if(y =='Disease') | (y == 'BiologicalProcess'):
                for z in hp_symptom_dict[key]["names"]:
                    if((y == 'Disease') & (len(UMLS) > 0)): 
                        a = ht.query(UMLS)[y]
                    else:
                        a = ht.query(z)[y]
                    # print(a)
                    for b in a: 
                        if b['name'].lower() == z:
                            # print('match')
                            # print(b)
                            # print(z)
                            try: 
                                fc = FindConnection(input_obj=b, output_obj='Gene', intermediate_nodes=None)
                                fc.connect(verbose=False)
                                df = fc.display_table_view()
                                # print("BD")
                                # print(df.shape[0])
                                if(df.shape[0] > 0):
                                    df = df[df["output_name"] != disease_name]
                                    edges_out_count = edges_out_count + df.shape[0]

                            except: 
                                print("Nope")
                            try: 
                                fc = FindConnection(input_obj=b, output_obj='Disease', intermediate_nodes=None)
                                fc.connect(verbose=False)
                                df = fc.display_table_view()
                                # print("BD")
                                # print(df.shape[0])
                                if(df.shape[0] > 0):
                                    df = df[df["output_name"] != disease_name]
                                    edges_out_count = edges_out_count + df.shape[0]
                            except: 
                                print("Nope")
        print("edges out")
        print(edges_out_count)
        hp_symptom_dict[key]["edges_out_count"] = edges_out_count
    return(hp_symptom_dict)


def get_similar_phen_indices(list1,list2,similarity,HP_dict):
    res = [] 
    i = 0
    while (i < len(list1)):
        append_i = False
        lookup = list1[i].lower()
        if('HP:' in list1[i]):
            if(list1[i]  in HP_dict):
                lookup = HP_dict[list1[i]]
        for j in list2:
                if(difflib.SequenceMatcher(None,lookup,j).ratio() > similarity):
                    append_i = True
        if(append_i): 
            res.append(i) 
        i += 1
    return(res)

def get_similar_bp_indices(list1,list2,similarity,go_dict):
    res = [] 
    i = 0
    while (i < len(list1)):
        append_i = False
        lookup = list1[i].lower()
        if('go:' in list1[i]):
            if list1[i] in go_dict:
                lookup = go_dict[list1[i]]
        for j in list2:
                if(difflib.SequenceMatcher(None,lookup,j).ratio() > similarity):
                    append_i = True
        if(append_i): 
            res.append(i) 
        i += 1
    return(res)

def get_similar_disease_indices(list1,list2,similarity):
    res = [] 
    i = 0
    while (i < len(list1)):
        append_i = False
        lookup = list1[i].lower()
        for j in list2:
                if(difflib.SequenceMatcher(None,lookup,j).ratio() > similarity):
                    append_i = True
        if(append_i): 
            res.append(i) 
        i += 1
    print(len(res))
    return(res)


def determined_genes_to_symptoms(gene_list, symptom_list):
    # gene -> phenotypic feature nodes
    print("Genes -> PhenotypicFeatures")

    df_list = []
    for x in gene_list: 
    #     print(x)
        try: 
            gene = ht.query(x)["Gene"][0]
            fc = FindConnection(input_obj=gene, output_obj='PhenotypicFeature', intermediate_nodes=None)
            fc.connect(verbose=False)
            df = fc.display_table_view()
            rows = df.shape[0]
            if(rows > 0):
                df_list.append(df)
        except:
            print(str(x) + " FAILED")
    if(len(df_list) > 0):
        top_gene_to_phenotypicFeature = pd.concat(df_list)

    ## Get names for HP ids
    HP_ids = top_gene_to_phenotypicFeature[top_gene_to_phenotypicFeature["output_name"].str.contains("HP:",regex=False)]["output_name"]
    HP_ids = list(HP_ids)
    HP_ids = list(dict.fromkeys(HP_ids))
    # len(HP_ids)
    HP_dict = {}
    for x in HP_ids: 
        HP_ID = x.split(':')[1]
        r = requests.get('https://biothings.ncats.io/hpo/phenotype/HP%3A' + HP_ID)
        res = r.json()
        if(('_id' in res) & ('name' in res)):
            HP_dict[res['_id']] = res['name'].lower()

    phen_indices = get_similar_phen_indices(list(top_gene_to_phenotypicFeature["output_name"]),symptom_list,0.95, HP_dict)
    phen_top = top_gene_to_phenotypicFeature.iloc[phen_indices,:]
    
    # phen_top = top_gene_to_phenotypicFeature
    # phen_top
    for index in range(phen_top.shape[0]):
    #     if("HP:" in row['output_name']):
    #     print(index)
        if(phen_top.iloc[index]["output_name"] in HP_dict):
            phen_top.iloc[index]["output_name"] = HP_dict[phen_top.iloc[index]["output_name"]]

    phen_top

    # gene -> bioprocess
    print("Genes -> Bioprocesses")
    df_list = []
    for x in gene_list: 
    #     print(x)
        try: 
            gene = ht.query(x)["Gene"][0]
            fc = FindConnection(input_obj=gene, output_obj='BiologicalProcess', intermediate_nodes=None)
            fc.connect(verbose=False)
            df = fc.display_table_view()
            rows = df.shape[0]
            if(rows > 0):
                df_list.append(df)
        except:
            print(str(x) + " FAILED")
    if(len(df_list) > 0):
        top_gene_to_bioprocesses = pd.concat(df_list)

    go_ids = top_gene_to_bioprocesses[top_gene_to_bioprocesses["output_name"].str.contains("go:",regex=False)]["output_name"]
    go_ids = list(go_ids)
    go_ids = list(dict.fromkeys(go_ids))
    # len(go_ids)
    go_dict = {}
    for x in go_ids: 
        go_ID = x.split(':')[1]
        r = requests.get('https://biothings.ncats.io/go_bp/geneset/GO%3A' + go_ID)
        res = r.json()
        if('name' in res):
            go_dict[res['_id']] = res['name'].lower()

    bp_indices = get_similar_bp_indices(list(top_gene_to_bioprocesses["output_name"]),symptom_list,0.95,go_dict)
    bioprocess_top = top_gene_to_bioprocesses.iloc[bp_indices,:]

    # bioprocess_top = top_gene_to_bioprocesses
    # Genes -> disease type "symptoms"
    print("Genes -> Diseases")
    df_list = []
    for x in gene_list: 
        try: 
            gene = ht.query(x)["Gene"][0]
            fc = FindConnection(input_obj=gene, output_obj='Disease', intermediate_nodes=None)
            fc.connect(verbose=False)
            df = fc.display_table_view()
            rows = df.shape[0]
            if(rows > 0):
                df_list.append(df)
        except:
            print(str(x) + " FAILED")
    if(len(df_list) > 0):
        top_gene_to_diseases = pd.concat(df_list)

    disease_indices = get_similar_disease_indices(list(top_gene_to_diseases["output_name"]),symptom_list,0.95)
    relevant_top_gene_to_diseases = top_gene_to_diseases.iloc[disease_indices,:]

    # relevant_top_gene_to_diseases = top_gene_to_diseases

    ## make dataframe with all genes -> symptoms
    all_gene_connections = pd.concat([bioprocess_top,phen_top,relevant_top_gene_to_diseases])
    all_gene_connections["output_name"] = all_gene_connections["output_name"].str.lower()
    return(all_gene_connections)

def get_gene_to_symptom_publication_counts(all_gene_connections):
    # get pulication counts for gene -> symptoms
    gene_to_symptom_pub_counts = {}
    for index, row in all_gene_connections.iterrows():
    #     if row["input_name"] in top_related_genes_covid_2_all_nodes_2_genes:
        current_pubcount = 0
        if(row["pred1_pubmed"] != None):
            current_pubcount = current_pubcount + row["pred1_pubmed"].count(",") + 1
        if row["input"] in gene_to_symptom_pub_counts:
            gene_to_symptom_pub_counts[row["input"]] = gene_to_symptom_pub_counts[row["input"]] + current_pubcount
        else: 
            gene_to_symptom_pub_counts[row["input"]] = current_pubcount
    return(gene_to_symptom_pub_counts)


def create_causes_dict(all_gene_connections):
    causes_df = all_gene_connections[all_gene_connections["pred1"] == "causes"]
    i = list(causes_df["input"])
    causes_dict = {x:i.count(x) for x in i}
    return(causes_dict)

# function that gets all connections to any node type from single gene node
def get_connection_normalizing_count(gene_list, node_type_list):
    # dictionary that keeps track of all connections from a gene to any node type 
    connection_dict = {}
    for gene_symbol in gene_list:
        gene_found = False
        gene_query = ht.query(gene_symbol)['Gene']
        for i in gene_query:
            if(i['SYMBOL'].lower() == gene_symbol.lower()):
                gene = i
                gene_found = True
        if(gene_found == True):
            count = 0
            input_object = gene
            for x in node_type_list:
                try: 
                    ## only look at direct connections
                    fc = FindConnection(input_obj=input_object, output_obj=x, intermediate_nodes=None)
                    fc.connect(verbose=False)
                    df = fc.display_table_view()
                    rows = df.shape[0]
                    count = count + rows
                except: 
                    print("gene " + str(gene_symbol) + " for node intermediate " + str(x) + " failed")
            connection_dict[gene_symbol]  = count
        else:
            print(gene_symbol + ' could not be found')
            connection_dict[gene_symbol] = 'Unknown'
    return(connection_dict)

def assemble_final_data_frame(all_gene_connections, connection_dict, sorted_disease_to_genes, sorted_disease_to_all_nodes_to_genes, top_two_step_genes_pub_counts, top_symptom_pub_counts, causes_dict, disease_symptoms_df):

    # make dictionary for final assembly of results
    results_dict = {}
    for i in range(all_gene_connections.shape[0]):
        if(all_gene_connections.iloc[i]["input"] in results_dict):
            results_dict[all_gene_connections.iloc[i]["input"]]["symptoms_associated"].append(all_gene_connections.iloc[i]["output_name"])
        else:
            results_dict[all_gene_connections.iloc[i]["input"]] = {
                "two_step_associations_to_genes" : sorted_disease_to_all_nodes_to_genes[all_gene_connections.iloc[i]["input"]] if all_gene_connections.iloc[i]["input"] in sorted_disease_to_all_nodes_to_genes else 0,
                "direct_associations_to_genes" : sorted_disease_to_genes[all_gene_connections.iloc[i]["input"]] if all_gene_connections.iloc[i]["input"] in sorted_disease_to_genes else 0,
                "symptoms_associated" : [all_gene_connections.iloc[i]["output_name"]]
            }



        
    # assmeble results into final dataframe
    dataframe_input = []
    for key in results_dict:

        ### Weight --- number - min / range


        connections_count = math.sqrt(connection_dict[key])
        # calculate "relevance_score" based on occurrences, publication counts, gene_normalizing counts 
        # relevance_score = ((((results_dict[key]["direct_associations_to_genes"]*10 
        #                     + results_dict[key]["two_step_associations_to_genes"]) 
        #                     * len(results_dict[key]["symptoms_associated"])*3) 
        #                     # + round(top_two_step_genes_pub_counts[key] / 5) 
        #                     # + round(top_symptom_pub_counts[key] / 5)
        #                     + (causes_dict[key] if key in causes_dict else 0)*20)
        #                     /connections_count)
        # assemble each row                                               
        current_result = {'gene': key,
                        "direct_disease_assoc": results_dict[key]["direct_associations_to_genes"], 
                        "two_step_assoc_to_disease": results_dict[key]["two_step_associations_to_genes"],
                        # "two_step_pub_count": top_two_step_genes_pub_counts[key] if key in top_two_step_genes_pub_counts else 0,
                        "disease_symptoms_gene_is_associated_with": results_dict[key]["symptoms_associated"],
                        "symptoms_associated_count": len(results_dict[key]["symptoms_associated"]),
                        # "disease_symptom_gene_pub_count": top_symptom_pub_counts[key],
                        # "causes_symptom_count": causes_dict[key] if key in causes_dict else 0,
                        "gene_connections_count": connection_dict[key]
                        # "relevance_score": relevance_score
                        }
        dataframe_input.append(current_result)
    
    final_df = pd.DataFrame(dataframe_input)
    
    
    symptom_score_dict = {}
    for index, row in disease_symptoms_df.iterrows():
        # print(row)
        for x in row["names"]:
            symptom_score_dict[x.lower()] = row["ISS"]

    # print(symptom_score_dict)
    final_symptom_scores = []

    for index, row in final_df.iterrows():
        current_score = 0
        current_symptoms = row["disease_symptoms_gene_is_associated_with"]
        for x in current_symptoms:
            singular = x.lower()[0:-1]
            plural = x.lower() + 's'
            if(x.lower() in symptom_score_dict):
                current_score = current_score + symptom_score_dict[x.lower()]
            elif(plural in symptom_score_dict):
                current_score = current_score + symptom_score_dict[plural]
            elif(singular in symptom_score_dict):
                current_score = current_score + symptom_score_dict[singular]
        final_symptom_scores.append(current_score)


    final_df["final_symptom_score"] = final_symptom_scores


    range_direct = max(list(final_df["direct_disease_assoc"])) - min(list(final_df["direct_disease_assoc"]))
    min_direct = min(list(final_df["direct_disease_assoc"]))

    range_two_step = max(list(final_df["two_step_assoc_to_disease"])) - min(list(final_df["two_step_assoc_to_disease"]))
    min_two_step = min(list(final_df["two_step_assoc_to_disease"]))

    range_symptoms = max(list(final_df["final_symptom_score"])) - min(list(final_df["final_symptom_score"]))
    min_symptoms = min(list(final_df["final_symptom_score"]))
    relevance_score = []
    for index, row in final_df.iterrows():
        current_direct = final_df["direct_disease_assoc"][index]
        current_two_step = final_df["two_step_assoc_to_disease"][index]
        current_symptom = final_df["final_symptom_score"][index]
        direct_weighted = (current_direct - min_direct)/range_direct
        two_step_weighted = (current_two_step - min_two_step)/range_two_step
        symptom_weighted = (current_symptom - min_symptoms)/range_symptoms
        
        current_rs = ((direct_weighted + two_step_weighted)*symptom_weighted)/math.sqrt(final_df["gene_connections_count"][index])
        relevance_score.append(current_rs)
    
    relevance_score_norm = [(float(i) - min(relevance_score))/(max(relevance_score)-min(relevance_score)) for i in relevance_score]
    final_df["relevance_score"] = relevance_score_norm
    # sort by relevance score
    final_df = final_df.sort_values(by=['relevance_score'], ascending=False)
    return(final_df)