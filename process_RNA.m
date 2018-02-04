close all
clc
clear
%%
d5fu_non=importdata('5fu_non.csv');
%%
d5fu_res=importdata('5fu_res.csv');
%%
drugname='5fu';
%%
xx=[d5fu_non,d5fu_res]';
x=xx;

%%
genenamed=importdata('gene_name.txt');
genename=genenamed';
marker=genename;
%%
name=importdata('pn.txt');
%%
y=importdata('y.txt');
%%
train_idx1 = randperm(25,19);
train_idx2 = randperm(35,26)+25;
train_idx=[train_idx1,train_idx2];
%%
test_idx=1:60;
test_idx(train_idx)=[];
%%