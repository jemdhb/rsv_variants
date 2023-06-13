#graphing and ds
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
#colors for graphs
c1="mediumvioletred"
c2=(0.990249903883122, 0.43280276816609, 0.20096885813148788)
c_pal=[c1,c2]
big_c_pal=[c1,"red",c2,"goldenrod"]

def set_seaborn_formatting(x=11.7,y=8.27):
    """clean up seaborn graphs
    """
    sns.set(rc={'figure.figsize':(x,y)})
    sns.set_theme(style="white")

def set_fig_labels(xlabel="",ylabel="",title="",fig=None):
    """set labels for figs pre/post graphing
    """
    plt.gca().set(xlabel=xlabel, ylabel=ylabel, title=title)
    #if fig has already been created
    if fig !=None:
        return_figure(fig)

def return_figure(plot):
    """Display clean version of graph for jupyter notebook
    """
    clean_names={"collection_month_year":"Collection date",
                 "ocurrence_counter":"Count",
                 "collection_year":"Collection year",
                 "count":"Count","Density":"Density",
                 "collection_date":"Collection date"}
    try:
        plot.set(xlabel=clean_names[plot.xaxis.get_label().get_text()],
                ylabel=clean_names[plot.yaxis.get_label().get_text()])
    except:
        pass

    fig=plot.get_figure()
    sns.despine(fig,top=True,right=True)
    plt.show()
    plt.close("all")
    return fig

def graph_barplot(df,x,y,pal):
    """Create custom barplot"""
    return return_figure(sns.barplot(data=df,x=x,
                                     y=y,palette=pal))

def graph_running_average_line(df,x, fig, color="blue",linewidth=2):
    """create running average to overlay onto another figure

    Args:
        df (dataframe): df storing info
        x (str): col name used for running avg
        fig (matplotlib figure): figure to overlay line on
        color (str, optional): determines line color. Defaults to "blue".
        linewidth (int, optional): determines line width. Defaults to 2.
    """
    value_counts_df=df[x].value_counts().rename_axis('unique_values').to_frame('counts').reset_index()
    value_counts_df.sort_values(["unique_values"],ascending=True,inplace=True)
    avg=[]
    n=0
    for item in value_counts_df["counts"]:
        if n==0:
            avg.append(item)
        else:
            new_avg = (avg[-1] * (n-1)/n)+ (item /n)
            avg.append(new_avg)
        n+=1
    sns.lineplot(x=value_counts_df["unique_values"],y=avg,color=color,
             alpha=0.7, linewidth=linewidth,linestyle='--')
    
def graph_countplot(df,x,average_line=True,bar_color="red",
                    line_color=None):
    """Create Countplot from df column x, w/ or w/out a running average line

    Args:
        df (df): Dataframe with our info
        x (str): df column name
        average_line (bool): Include average line on countplot
        bar_color (str, optional): Colors countplot. Defaults to "red".
        line_color (str, optional): Colors average line. Defaults to None.
    """
    order=None
    if isinstance(df, pd.DataFrame)==True:
        order=sorted(set(df[x]))
        ax=sns.countplot(data=df,x=x,order=order,color=bar_color,alpha=0.3)
    #if graphing a list, not a df col
    else:
        ax=sns.countplot(x=x,order=order,color=bar_color,alpha=0.3)

    if average_line:
        graph_running_average_line(df,x,ax,color=line_color,linewidth=4)
    return return_figure(ax)

def graph_lineplot(df,x,y,hue,style,pal,linewidth=3,alpha=0.7,linestyle="-",
                   order=["B","A"]):
    """create custom lineplot that splits lines and coloring by order (variant
    name)
    """
    ax=sns.lineplot(data=df, x=x,y=y,hue=hue,style=style,palette=pal,
                 linewidth=linewidth,alpha=alpha,linestyle=linestyle,
                 style_order=order)
    return return_figure(ax)

def graph_kde(df,x,hue,pal,start,end,alpha=0.7,multiple="fill"):
    """graph a custom kde with the midpoint annotated
    """
    ax=sns.kdeplot(data=df, x=x, hue=hue, palette=pal,alpha=alpha,
                multiple=multiple)
    #halfway point line
    plt.axhline(y=0.5, color='black', linestyle='--',alpha=0.7)
    #trim malformed dates
    plt.xlim(start,end)
    return return_figure(ax)

def graph_scatterplot(df,x,y,alpha=0.3,color="black",pal=None,hue=None,
                      size=None):
    """graph custom scatterplot
    """
    ax=sns.scatterplot(df,x=x,y=y,alpha=alpha,color=color,palette=pal,hue=hue,
                       size=size, legend=size==None)
    return return_figure(ax)

def graph_stripplot(df,x,color,hue=None,palette=None,alpha=0.3):
    """create stripplot

    Args:
        df (df): dataframe with our info to graph
        x (str): df column with relevant data
        color: default point colors
        hue (str, optional): column to differentiate points. 
        palette (str, optional): color palette to differentiate points. 
        alpha (float, optional): Determines opacity of dots 0-1.
    """
    #if graphing a list
    if isinstance(df, pd.DataFrame)==False:
        return_figure(sns.stripplot(data=None,
                                    x=x,size=2,alpha=alpha,
                                    color=color,hue=hue,palette=palette))
    #if graphing a df
    else:
        df[x]=pd.to_datetime(df[x])
        return_figure(sns.stripplot(data=df, x=x,size=2,alpha=alpha,
                                    color=color,hue=hue,palette=palette))

