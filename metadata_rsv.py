from datetime import datetime
import pandas as pd
import numpy as np
from collections import Counter
from statistics import mean
#for embedded bash commands
import subprocess
#classification and regression
from sklearn.model_selection import GridSearchCV, train_test_split
from sklearn.linear_model import LogisticRegression, LinearRegression
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import accuracy_score,precision_score, recall_score, f1_score
import warnings
#ordering for metadata extraction
variant_name_i=1
location_name_i=2
iep=3

def get_names_simple(l):
    """shorten variant names to (RSV)A or (RSV)B
    """
    return [item[find_nth(item, "/", 1)+1:find_nth(item, "/", 2)]
            for item in l]

def country_to_continent(df):
    keys=open("blast_results/countries.txt","r").read().split("\n")
    values=open("blast_results/continents.txt","r").read().split("\n")
    res = dict(zip(keys, values))
    df["Continent"]=[res[item] for item in df["Country"]]
    continent_df=df.groupby("Continent").sum().sort_values(by="Submission count",
                                                           ascending=False)
    return continent_df

def get_categorical_mean(df,col):
    """Get the mean value for col in df
    """
    categorical_counts=list(df[col].value_counts())
    mean_categorical_counts=mean(categorical_counts)
    return mean_categorical_counts

def extract_metadata():
    """grab headers from fasta
    """
    command = """cat gisaid_rsv_2023_04_13_05.fasta | grep ">" | cut -d ">" -f2- >metadata.txt"""
    process = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()

def txt_to_csv(file:str, delim:str):
    """make fasta metadata pandas readable
    """
    #create csv with header
    output=file[:file.rfind(".")]+".csv"
    output_file=open(output,"w")
    output_file.write("variant_name,epi_isl,collection_date\n")
    input_file=open(file,"r")
    for line in input_file:
        #only care about headers in txt file
        if delim in line:
            output_file.write(line.replace(delim,",").strip()+"\n")  
    
    input_file.close()
    output_file.close()

    return output

def find_nth(string, substring, n):
   """used in extract_fasta_header"""
   if n==0:
       return 0
   if (n == 1):
       return string.find(substring)
   else:
       return string.find(substring, find_nth(string, substring, n - 1) + 1)

def extract_fasta_header(df,col,n, delim="/"):
    """Extract header info stored in df col, trimming to delim n"""
    start=n
    end=n+1
    last_delim=3
    if n>last_delim:
        #use rfind for ending delims
        return [item[item.rfind(delim)+1:] for item in df[col]]
    #otherwise recurse
    return [item[find_nth(item, delim, start)+1:find_nth(item, delim, end)] for\
            item in df[col]]
    
def variants_over_time(input_file:str):
    """rsv information driver

    Args:
        input_file (str): metadata file
    """
    #read in metadata
    rsv_df=pd.read_csv(txt_to_csv(input_file,"|"))
    
    subtype_extraction=extract_fasta_header(rsv_df,"variant_name",variant_name_i)
    
    #counts of every country that submitted RSV data
    rsv_by_country_list=extract_fasta_header(rsv_df,"variant_name",
                                                location_name_i)
    #clean up groupings
    rsv_by_country_list=[item.replace("_"," ").replace(
                        "England","United Kingdom").replace(
                        "NewZealand","New Zealand")
                        for item in rsv_by_country_list]
    rsv_by_country=Counter(rsv_by_country_list)
    #force rsv_by_country in df form
    rsv_by_country_df=pd.DataFrame.from_dict(
        rsv_by_country,orient="index").reset_index().sort_values(
        by=0,ascending=False)
    rsv_by_country_df=rsv_by_country_df.rename(columns=
        {"index":"Country",0:"Submission count"})
    
    isp_num=Counter(extract_fasta_header(rsv_df,"variant_name",iep))

    #counts of every year with RSV data
    collection_year=Counter(extract_fasta_header(rsv_df,"variant_name",4))
    #force collection_year into df form
    collection_year_df=pd.DataFrame.from_dict(
        collection_year,orient="index").reset_index().sort_values(
        by="index",ascending=True)
    
    return rsv_df,subtype_extraction,rsv_by_country_df,isp_num,collection_year_df

def grid_search(classifier, l,X_train,y_train):
    """generalized classification optimization

    Args:
        classifier (classification constructor): 
        l (str): generic params to optimize on
        X_train (df): data to base classification on 
        y_train (df): known classification results

    Returns:
        _type_: _description_
    """
    if classifier==LinearDiscriminantAnalysis():
        search = GridSearchCV(classifier, l, scoring="accuracy")

    else:
        search = GridSearchCV(classifier, l, cv=3)

    search.fit(X_train, y_train)
    res=[]
    for item in l.keys():
        res.append(search.best_params_.get(item))
        l[item]=search.best_params_.get(item)
    return res

def calculate_scores(clf,X_train, y_train,X_test,y_test):
    """classify our data with the inputted classifier and then rate the 
    classifiers performance
    """
    #fit classifier to training data
    clf.fit(X_train, y_train)
    #for future results presentation
    name = clf.__class__.__name__
    train_predictions = clf.predict(X_test)
    # calculate accuracy, precision, recall, and fscore keeping first 4 decimals
    acc = ("{:.2%}".format(accuracy_score(y_test, train_predictions)))
    precision = ("{:.2%}".format(precision_score(y_test, train_predictions, average = 'macro'))) 
    recall =  ("{:.2%}".format(recall_score(y_test, train_predictions, average = 'macro')))
    f_score = ("{:.2%}".format(f1_score(y_test, train_predictions, average = 'macro')))
    return [name,acc,precision,recall,f_score]


        
if __name__ == "__main__":
    variants_over_time('metadata.txt')