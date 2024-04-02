#!/usr/bin/env python 
# -*- coding: utf-8 -*- 
#Usage: python dotplot.py AHC.gff AAS.gff AHC.lens AAS.lens GD_pairs.txt
#Usage: 20240310, Yiyogn Zhao modified

#AHC.gff
#for gff, the last column is the gene count for each chromosome.
#chr1	AHC_0004600	12837	26777	+	1
#chr1	AHC_0004073	33171	35791	+	2
#chr1	AHC_0004946	46794	47258	-	3

#AHC.lens
#chr1	23021215
#chr2	18778401
#chr3	19321588

#AAS.gff
#for gff, the last column is the gene count for each chromosome.
#Chr1	AAS_0012495	3631	5899	+	1
#Chr1	AAS_0011322	5928	8737	-	2
#Chr1	AAS_0007791	11649	13714	-	3

#AAS.lens
#Chr1	30425192
#Chr2	19696821
#Chr3	23459804

#GD_pairs.txt
#if first iput is AHC.gff, so the GD pairs first column is AHC genes; 
#AHC_0000001	AAS_0011214	NONE
#AHC_0000001	AAS_0011214	NONE
#AHC_0000001	AAS_0011214	NONE

import matplotlib 
matplotlib.use('Agg') 
import sys 
import re 
import numpy as np 
import matplotlib.pyplot as plt; plt.rcdefaults() 
import matplotlib.patches as mpatches 
import time 
import gc


def read_gff(fn):
    f=open(fn) 
    data,dict=[],{} 
    for line in f.readlines(): 
        a=line.strip().split("\t") 
        dict[a[1]]=[a[0],a[2],a[3],a[4],a[5]]
        data.append(a) 
    return data,dict
def read_lens(fn): 
    fp=open(fn) 
    data=[] 
    for row in fp.readlines(): 
        r1,r2=row.split() 
        data.append([r1,r2]) 
    return data
def plot_chr1(lens,gl,gl2,mark,name): 
    total_lens=sum([float(k[1]) for k in lens]) 
    step=gl/float(total_lens) 
    gl_start,n,start_x=0.95,0,0.05 
    mark_y=0.04 
    align = dict(family='Times New Roman',style='normal',horizontalalignment ="center", verticalalignment="center") 
    for k in lens: 
        n+=float(k[1]) 
        mark_new=str(mark)+str(k[0]) 
        x=gl_start-float(n)*step 
        mark_x=x+0.5*float(k[1])*step 
        plt.plot([start_x,start_x+gl2],[x,x],linestyle = '-',color='black', linewidth=0.5) 
        plt.text(mark_y,mark_x,mark_new,color='black',fontsize = 12,rotation = 90,weight='semibold',**align) 
    plt.plot([start_x,start_x+gl2],[gl_start,gl_start],linestyle='-', color='black',linewidth=1) 
    plt.text(mark_y-0.02,0.5*(2*gl_start-gl),name,color='black',fontsize = 18,rotation = 90,weight='semibold',**align) 
    return step

def plot_chr2(lens,gl,gl2,mark,name): 
    total_lens=sum([float(k[1]) for k in lens]) 
    step=gl/float(total_lens) 
    gl_start,n,start_x=0.05,0,0.95 
    mark_y=0.96 
    align = dict(family='Times New Roman',style='normal',horizontalalignment="center", verticalalignment="center") 
    for k in lens: 
        n+=float(k[1]) 
        mark_new=str(mark)+str(k[0]) 
        x=gl_start+float(n)*step 
        mark_x=x-0.5*float(k[1])*step 
        plt.plot([x,x],[start_x,start_x-gl2],linestyle='-',color='black', linewidth=0.5) 
        plt.text(mark_x,mark_y,mark_new,color='black',fontsize = 12,rotation = 0,weight='semibold',**align) 
    plt.plot([gl_start,gl_start],[start_x,start_x-gl2],linestyle='-', color='black',linewidth=1) 
    plt.text(0.5*(2*gl_start+gl),mark_y+0.02,name,color='black',fontsize = 18,rotation = 0,weight='semibold',**align) 
    return step
def gene_location(gff, lens, step): 
    loc_gene,dict_chr,n={},{},0 
    for i in lens: 
        dict_chr[i[0]]=n 
        n+=float(i[1]) 
    for k in gff: 
        if k[0] not in dict_chr.keys(): 
            continue 
        loc=(float(dict_chr[k[0]])+float(k[3]))*step 
        loc_gene[k[1]]=loc 
    return loc_gene

def read_gd_pairs(fn):
	f=open(fn)
	dict_gd,dict_gd1={},{}
	for line in f.readlines():
		a=line.strip().split("\t")
		dict_gd1[str(a[0])+":"+str(a[1])]=str(a[2])
	for pair, wgd_notch in dict_gd1.items():
		fd=pair.strip().split(":")
		gene1=fd[0]
		gene2=fd[1]
		index= gene1+":"+gene2
		reverse_index=gene2+":"+ gene1
		dict_gd[index]=dict_gd1[index]
		if reverse_index not in dict_gd1.keys():
			dict_gd[reverse_index]=dict_gd1[index]
	return dict_gd
    
def plot_dot(root,data,loc1,loc2,gl):
	gl_start1,gl_start2=0.95,0.05
	for pair, wgd_notch in dict_gd.items():
		fd=pair.strip().split(":")
		gene1=fd[0]
		gene2=fd[1]
		index= gene1+":"+gene2
		if gene1 in loc1.keys() and gene2 in loc2.keys():
			x,y=loc1[gene1],loc2[gene2]
			x,y=gl_start1-x,gl_start2+y
			if index in dict_gd.keys():
				if dict_gd[index] == "red":
					DrawCircle(root,[y,x],0.001,'red',0.6)
				if dict_gd[index] == "green":
					DrawCircle(root,[y,x],0.001,'green',0.6)
				if dict_gd[index] == "orange":
					DrawCircle(root,[y,x],0.001,'orange',0.6)
				if dict_gd[index] == "blue":
					DrawCircle(root,[y,x],0.001,'blue',0.6)
				if dict_gd[index] == "NONE":
					DrawCircle(root,[y,x],0.001,'gray',0.6)   

def DrawCircle(ax, loc, radius,color, alpha):
    circle = mpatches.Circle(loc, radius, edgecolor="none", facecolor=color, alpha=alpha)
    ax.add_patch(circle)

if __name__=="__main__":
    gff1=sys.argv[1]
    gff2=sys.argv[2]
    lens1=sys.argv[3]
    lens2=sys.argv[4]
    gd_pairs=sys.argv[5]
    spe1=sys.argv[6]
    spe2=sys.argv[7]
    plt.figure(figsize=(10, 10)) 
    root = plt.axes([0, 0, 1, 1]) 
    align = dict(family='Arial',style='normal',horizontalalignment="center", verticalalignment="center") 
    t1=time.time() 
    print(f"Dotplot {gd_pairs} are ready to begin") 
    gff_1,dict_gff1=read_gff(gff1) 
    gff_2,dict_gff2=read_gff(gff2) 
    t2=time.time() 
    print("Reading gff took "+str(t2-t1)+" second") 
    lens_1=read_lens(lens1)
    lens_2=read_lens(lens2)
    t3=time.time() 
    print("Reading lens took "+str(t3-t2)+" second") 
    gl1,gl2=0.92,0.92 
    step_1=plot_chr1(lens_1,gl1,gl2,'',spe1) 
    step_2=plot_chr2(lens_2,gl2,gl1,'',spe2) 
    dict_gd=read_gd_pairs(gd_pairs)
    t4=time.time() 
    print(f"Reading {gd_pairs}_pairs took "+str(t4-t3)+" second")  
    gene_loc_1=gene_location(gff_1,lens_1,step_1)
    gene_loc_2=gene_location(gff_2,lens_2,step_2) 
    t5=time.time() 
    print(f"Dealing {gd_pairs}_file took "+str(t5-t4)+" second") 
    gc.collect() 
    plot_dot(root,dict_gd,gene_loc_1,gene_loc_2,gl1) 
    t6=time.time() 
    print("Ploting dot took "+str(t6-t5)+" second") 
    root.set_xlim(0, 1) 
    root.set_ylim(0, 1) 
    root.set_axis_off() 
    plt.savefig(gd_pairs+"_dotplot.pdf",dpi=500) 
    plt.savefig(gd_pairs+"_dotplot.png",dpi=500) 
    t7=time.time() 
    print(f"Dotplot {gd_pairs} totaly took "+str(t7-t1)+" second")


