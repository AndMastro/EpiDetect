import networkx as nx
import sys
import csv

if __name__ == "__main__":
    

    arg_len = len(sys.argv)

    if arg_len < 3:
        print("ERROR: not enough arguments. Specify file and centrality measure")
        sys.exit()

    DATA_FILE = sys.argv[1]
    centrality = sys.argv[2]


    edges = []

    with open(DATA_FILE, "r") as dataFile:
        for row in dataFile:
            line = row.strip().split("\t")
            edges.append((line[2], line[3]))
    
    netName = DATA_FILE
    G = nx.Graph(name = DATA_FILE)
    G.add_edges_from(edges)

    print(nx.info(G))

    if centrality == "degree":

        avg_degree = sum(dict(G.degree()).values())/G.number_of_nodes()
        print("Average " + centrality + ": " + str(avg_degree))

        deg = G.degree()

        deg_dict = {}
        for gene in deg:
            deg_dict[gene[0]] = gene[1]

        # Create graph usong only central genes (deg > avg_degree)


        central_genes = []
        central_measures = []
        for pair in deg:
            if pair[1] > avg_degree:
                central_genes.append(pair[0])
                central_measures.append(pair[1])
        
        central_genes = [x for _,x in sorted(zip(central_measures,central_genes), reverse=True)]
        central_edges = []

        for gene1, gene2 in edges:
            if gene1 in central_genes and gene2 in central_genes:
                central_edges.append((gene1, gene2))

        central_G = nx.Graph(name = "Central Genes Network")

        central_G.add_edges_from(central_edges)

        print(nx.info(central_G))
        

        SAVE_PATH = DATA_FILE[:-4] + "_" + centrality +"_CENTRAL_GENES.txt"
        
        with open(SAVE_PATH, "w+") as saveFile:
            saveFile.write("Avg " + centrality +": " + str(avg_degree) + "\n\n")
            saveFile.write("Gene\tDegree\n")
            for gene in central_genes:
                saveFile.write(gene + "\t" + str(deg_dict[gene]) + "\n")
        

        '''deg_list = sorted(list(G.degree()), reverse = True, key = lambda x:x[1])

        central_deg_list = []

        for gene in deg_list:
            if gene[0] in genes_list:
                central_deg_list.append(gene)
                
        print(sorted(central_deg_list, reverse = True, key = lambda x:x[1]))'''
        
    else:
        centrality_measure = None
        avg_centrality = None

        if centrality == "betweenness":

            centrality_measure = nx.betweenness_centrality(G)
            avg_centrality = sum(centrality_measure.values())/G.number_of_nodes()

        elif centrality == "closeness":

            centrality_measure = nx.closeness_centrality(G)
            avg_centrality = sum(centrality_measure.values())/G.number_of_nodes()

        elif centrality == "eigenvector":
            centrality_measure = nx.eigenvector_centrality(G)
            avg_centrality = sum(centrality_measure.values())/G.number_of_nodes()

        elif centrality == "eccentricity":
            centrality_measure = nx.eccentricity(G)
            avg_centrality = sum(centrality_measure.values())/G.number_of_nodes()

        elif centrality == "diameter":
            diameter = nx.diameter(G)
            print("\nNetwork diameter (max eccentricity): " + str(diameter))
            sys.exit()

        elif centrality == "average_distance":
            distance = nx.average_shortest_path_length(G)
            print("\nNetwork average distance (avg shortest path): " + str(distance))
            sys.exit()
        
        elif centrality == "radiality":
            print("\nNOT YET IMPLEMENTED")
            sys.exit()
        
        elif centrality == "stress":
            print("\nNOT YET IMPLEMENTED")
            sys.exit()

        elif centrality == "centroid_value":
            print("\nNOT YET IMPLEMENTED")
            sys.exit()

        elif centrality == "bridging":
            print("\nNOT YET IMPLEMENTED")
            sys.exit()

        elif centrality == "edge_betweenness":
            print("\nNOT YET IMPLEMENTED")
            sys.exit()
            
        else:
            print("Centrality measure not implemented: choose from [degree, betweenness, closeness, eigenvector, eccentricity, diameter, average_distance, radiality, stress, centroid_value, bridging, edge_betweenness]")
            sys.exit()

        print("Average " + centrality + ": " + str(avg_centrality))

        # Create graph usong only central genes (deg > avg_degree)


        central_genes = []
        central_measures = []
        for gene in centrality_measure:
            if centrality_measure[gene] > avg_centrality:
                central_genes.append(gene)
                central_measures.append(centrality_measure[gene])
        
        central_genes = [x for _,x in sorted(zip(central_measures,central_genes), reverse=True)]
        central_edges = []
        
        for gene1, gene2 in edges:
            if gene1 in central_genes and gene2 in central_genes:
                central_edges.append((gene1, gene2))
        central_G = nx.Graph(name = "Central Genes Network")

        central_G.add_edges_from(central_edges)

        print(nx.info(central_G))

        central_G.nodes()

        SAVE_PATH = DATA_FILE[:-4] + "_" + centrality +"_CENTRAL_GENES.txt"
        

        with open(SAVE_PATH, "w+") as saveFile:
            saveFile.write("Avg " + centrality +": " + str(avg_centrality) + "\n\n")
            saveFile.write("Gene\t"+ centrality +"\n")
            for gene in central_genes:
                saveFile.write(gene + "\t" + str(centrality_measure[gene]) + "\n")

        
    # To add: take in input interaction file and centrality measure to use. Then, return a file with central genes: GENE_NAME DEG_ORIGINAL DEG_CENTRAL_NETWORK + avg degree in both networks (may be useful). To that for now, then, decide if consider edge weight.
