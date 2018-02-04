% mas5 data
close all
clc
clear
warning('off','all')
%%
load drug_RNA_5fu.mat
%%
a=(~x);
b=sum(a);
x=x(:,b<45);
genename=genename(b<45);
marker=marker(b<45);
y_drug = y; good_idx = ~isnan(y_drug);
feature_left_fix=100;
%%
name = name(good_idx);
x = x(good_idx,:);
y_drug = y_drug(good_idx);

%%
x_select = x;
y_select = y_drug;
name_select = name;

pos_z = 26:60; 
neg_z = 1:25;

%%
remain = 1:size(x,2);
sorted_feature = [];
n_feature = 0; step = 100;
label_select=logical(y);
%%
while n_feature < size(x,2)
    svm = svmtrain(x_select(train_idx,remain), label_select(train_idx),'kernel_function','linear','boxconstraint',1,'method','SMO');
    w = svm.Alpha'*svm.SupportVectors;
    if (length(remain) > step) 
        [~,I]=sort(abs(w)); I = I(step:-1:1);
    else
        [~,I]=min(abs(w)); I = I(1);
    end
    
    next_feature = remain(I);
    sorted_feature = [next_feature sorted_feature];
    remain(I) = [];
    
    n_feature = length(sorted_feature);
    if (n_feature == size(x,2)) || (mod(n_feature,step) == 0)
        fprintf('%d ', n_feature);
        if mod(n_feature/step,10) == 0, fprintf('\n'); end
    end
end
fprintf('\n');

%% 
    %test the features left on testing data
    acc = [];
    for f = [1:step]
        xtrain = x_select(train_idx,sorted_feature(1:f)); 
        ytrain = label_select(train_idx);
        xtest = x_select(test_idx,sorted_feature(1:f)); 
        ytest = label_select(test_idx);
        cp = classperf(label_select,'Positive',1,'Negative',0);
        svm = svmtrain(xtrain, ytrain,'kernel_function','linear','boxconstraint',1,'method','SMO');
        classOut = svmclassify(svm,xtest);
        cp = classperf(cp,classOut,test_idx);
        acc = [acc [cp.CorrectRate cp.Sensitivity cp.Specificity]'];

        if (f == size(x,2)) || (mod(f,1) == 0)
            fprintf('%d ', f);
            if mod(f/1,10) == 0, fprintf('\n'); end
        end
    end
    fprintf('\n');
%%  
    %make sure number seleted features is greater than threshold
    [~,I] = max(acc(1,:)); 
    if I(1)<feature_left_fix
        I=feature_left_fix;
    else
        I=I(1);
    end
    max_feature = sorted_feature(1:I); accuracy = acc(1,I)*100;
    feature_name = marker(max_feature);
    gene_name = genename(max_feature);
%%
    %Leave One Out
    lsa=[];
    cp = classperf(label_select,'Positive',1,'Negative',0);
    for i=1:size(x,1)
        loo_test_idx = false(size(x,1),1); loo_test_idx(i) = true;
        loo_train_idx = ~loo_test_idx;
        loo_xtrain = x(loo_train_idx,max_feature); 
        loo_ytrain=label_select(loo_train_idx);
        loo_xtest = x(loo_test_idx,max_feature);
        svm = svmtrain(loo_xtrain, loo_ytrain,'kernel_function','linear','boxconstraint',1,'method','SMO');
        classOut = svmclassify(svm,loo_xtest);
        cp = classperf(cp,classOut,loo_test_idx);
        loow = svm.Alpha'*svm.SupportVectors;
        loow_noscale = loow .* svm.ScaleData.scaleFactor;
        loob_noscale = sum(loow .* svm.ScaleData.scaleFactor .* svm.ScaleData.shift) + svm.Bias;
        looscore = loo_xtest*loow_noscale' + loob_noscale; looscore = -looscore;
        lsa=[lsa,looscore];
    end
    %%loo_acc = [cp.CorrectRate cp.Sensitivity cp.Specificity] * 100;
%%    
    %print out leave one out result
    X = x_select(:,max_feature); Y = label_select;
    xtrain = x_select(train_idx,max_feature); ytrain = label_select(train_idx);
    svm = svmtrain(xtrain, ytrain,'kernel_function','linear','boxconstraint',1,'method','SMO');
    w = svm.Alpha'*svm.SupportVectors;
    XX = bsxfun(@plus, X, svm.ScaleData.shift); XX = bsxfun(@times, XX, svm.ScaleData.scaleFactor);
    (XX*w'+svm.Bias)' < 0 == Y';
    w_noscale = w .* svm.ScaleData.scaleFactor;
    b_noscale = sum(w .* svm.ScaleData.scaleFactor .* svm.ScaleData.shift) + svm.Bias;
    (X*w_noscale'+b_noscale)' < 0 == Y';
    assert( all((X*w_noscale'+b_noscale < 0) == (XX*w'+svm.Bias < 0) ) )   
    score = x(:,max_feature)*w_noscale' + b_noscale; score = -score;
    pos_idx = 26:60; 
    neg_idx = 1:25; 
    %tmp = 1:sum(select_idx);
    %train_mark = false(size(x,1),1); 
    %train_mark(tmp(train_idx)) = true;
    %thry = 0.5*std_y;
    output_file = ['svm4_' drugname];
    %%
    score_t=lsa;
    p_label=score_t>0;
    tp=length(intersect(find(p_label>0),find(label_select>0)));
    tn=length(intersect(find(p_label<1),find(label_select<1)));
    fp=length(intersect(find(p_label>0),find(label_select<1)));
    fn=length(intersect(find(p_label<1),find(label_select>0)));
    Sensitivity=tp/(tp+fn)*100;
    Specificity=tn/(fp+tn)*100;
    label_select=label_select+0;
    p_label=p_label+0;
    test_acc=(1-nansum(abs(p_label-label_select'))/length(p_label))*100;
    %%
    pos_idx_l=find(label_select>0);
    neg_idx_l=find(label_select<1);
    clf, hold on
    h(1)=plot( lsa(pos_idx_l), y_drug(pos_idx ) , 'or' ,'MarkerSize',15,'LineWidth',3,'MarkerFaceColor', 'r');
    h(2)=plot( lsa(neg_idx_l), y_drug(neg_idx ) , 'ob' ,'MarkerSize',15,'LineWidth',3,'MarkerFaceColor', 'b');
    %psx=score_t(~train_mark );
    %psy=y_drug(~train_mark );
    %psl=label_select(~train_mark);
    %lst=length(psx);
    %for ki=1:lst
    %    if(psl(ki)>0)
    %        plot( psx(ki), psy(ki) , 'or' ,'MarkerSize',15,'LineWidth',3,'MarkerFaceColor', 'r');
    %    else
    %        plot( psx(ki), psy(ki) , 'ob' ,'MarkerSize',15,'LineWidth',3,'MarkerFaceColor', 'b');
    %    end
    %end
    
    %h(13)=plot( [min(lsa)-.3 max(lsa)+.3],[mean(y_drug)-thry, mean(y_drug)-thry],':', 'Color', [0.5 0.5 0.5],'LineWidth',3 );
    %h(15)=plot( [min(lsa)-.3 max(lsa)+.3],[mean(y_drug)+thry, mean(y_drug)+thry],':', 'Color', [0.5 0.5 0.5],'LineWidth',3 );
    set(h(1), 'MarkerFaceColor', 'r');
    set(h(2), 'MarkerFaceColor', 'b');
    t(1)=title(sprintf('%s\nacc=%.2f%% sn=%.2f%% sp=%.2f%%', drugname, [test_acc,Sensitivity,Specificity]));
    t(2)=xlabel( 'Predicted score');
    t(3)=ylabel( 'True label');
    set(t, 'FontSize',30);
    set(gca,'FontSize',30);
    xlim([min(lsa)-.3, max(lsa)+.3]);
    ylim([min(y_drug)-.3, max(y_drug)+.3]);
    legend('Sensitive','Resistant', 'Location','bestoutside');
    
    h_vert = line([0 0], ylim);
    set(h_vert, 'LineStyle', '--', 'Color', 'k', 'LineWidth', 2);
    h_horiz = line(xlim, [1 1]*mean(y_drug));
    set(h_horiz, 'LineStyle', '--', 'Color', 'k', 'LineWidth', 2);
    lx = xlim; ly = ylim;
    th(1) = text(lx(1), ly(1), 'True Negative', 'HorizontalAlignment', 'Left', 'VerticalAlignment', 'Bottom');
    th(2) = text(lx(1), ly(2), 'False Negative', 'HorizontalAlignment', 'Left', 'VerticalAlignment', 'Top');
    th(3) = text(lx(2), ly(2), 'True Positive', 'HorizontalAlignment', 'Right', 'VerticalAlignment', 'Top');
    th(4) = text(lx(2), ly(1), 'False Positive', 'HorizontalAlignment', 'Right', 'VerticalAlignment', 'Bottom');
    set(th, 'FontSize', 20);
    print_figure(gcf,[15 10], output_file,'-dpdf')