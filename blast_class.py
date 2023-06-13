import pandas as pd
import metadata_rsv as mr

def create_blast_cloud_from_homology_results(file_name,sep,one_node=False):
    """From homology results, create blast nodes which are then grouped into
    a blast cloud

    Args:
        file_name (str): file name of homology results
        sep (str): separator to split data
        one_node (bool, optional): If theres only one node dont make a list.
        Defaults to False.
    """
    file=open(file_name)
    chunk=""
    node_list=[]
    for line in file:
        if sep in line and chunk!="":
            node_list.append(feed_into_blast_node(chunk))
            chunk=""
            if one_node:
                return blast_cloud(node_list)
        chunk+=line
    return blast_cloud(node_list)

def feed_into_blast_node(blast_chunk):
    """Feed part of ncbi homology results to create a blast node

    Args:
        blast_chunk (str): str chunk with our relevant blast node info

    Returns:
        blast_node: blast node class which contains relevant results 
    """
    blast_chunk_list=blast_chunk.split("\n")
    #variant name
    name=blast_chunk_list[0]
    #alignment length
    length=blast_chunk_list[1]
    length=int(length[length.index("=")+1:].strip())
    #score and e value of alignment
    score_eval_list= blast_chunk_list[2].split(",")
    score,eval=score_eval_list[0],score_eval_list[1]
    score=int(score[score.index("(")+1:score.index(")")].strip())
    eval=float(eval[eval.index("=")+1:].strip())
    #homology and gap% of alignment
    identity_gap_list=blast_chunk_list[3].split(",")

    homology=identity_gap_list[0]
    homology=homology[homology.index("=")+1:homology.index("(")].strip()
    homology=int(homology[:homology.index("/")])/int(homology[homology.index("/")+1:])

    gap_perc=identity_gap_list[1]
    gap_perc=gap_perc[gap_perc.index("=")+1:gap_perc.index("(")].strip()
    gap_perc=int(gap_perc[:gap_perc.index("/")])/int(gap_perc[gap_perc.index("/")+1:])
    return blast_node(name,length,homology,gap_perc,eval,score)

class blast_node:
    """class to parse ncbi homology results
    """
    #node constructor: will only ever init when knowing all info
    def __init__(self,v_name,v_len,v_identity,v_gap_perc,v_eval,v_score):
        #name of variant
        self.name=v_name
        #sequence length of variant
        self.length=v_len
        #homology of variant to est. reference
        self.homology=v_identity
        #% of gaps in alignment to reference
        self.gap_perc=v_gap_perc
        #evalue of ncbi alignment
        self.eval=v_eval
        #alignment score of ncbi alignment
        self.score=v_score
    def get_name(self):
        return self.name
    def get_length(self):
        return self.length
    def get_homology(self):
        return self.homology
    def get_gap_perc(self):
        return self.gap_perc
    def get_eval(self):
        return self.eval
    def get_score(self):
        return self.score
    #print all relevant data
    def get_vals(self):
        print("name:",self.get_name())
        print("sequence length:",self.get_length())
        print("homology:",self.get_homology())
        print("gaps:",self.get_gap_perc())
        print("evalue",self.get_eval())
        print("score",self.get_score())
    
class blast_cloud:
    """group of blast nodes
    """
    def __init__(self,node_list):
        self.node_list=node_list  
    
    def get_cloud_names_simple(self):
        """get shortened variant names to more easilt group related variants
        """
        return mr.get_names_simple(self.get_cloud_names())
 
    def get_cloud_names(self):
        return [item.name for item in self.node_list] 
    
    def get_cloud_lengths(self):
        return [item.length for item in self.node_list]
    
    def get_cloud_homologies(self):
        return [item.homology for item in self.node_list]
    
    def get_cloud_gap_percs(self):
        return [item.gap_perc for item in self.node_list]
    
    def get_cloud_evals(self):
        return [item.eval for item in self.node_list]
    
    def get_cloud_scores(self):
        return [item.score for item in self.node_list]
    
    def get_cloud_df(self):
        """turn blast_cloud into df for graphing purposes
        """
        class_df=pd.DataFrame(
            {"homology":self.get_cloud_homologies(),
             "seq_len":self.get_cloud_lengths(),
             "gap_perc":self.get_cloud_gap_percs(),
             "evalue":self.get_cloud_evals(),
             "score":self.get_cloud_scores(),
             "name":self.get_cloud_names_simple()})
        return class_df