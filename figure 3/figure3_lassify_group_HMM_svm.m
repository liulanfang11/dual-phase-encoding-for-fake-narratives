% using the brain state features derived from HMM analysis to 
% decode participants' condition label.

clear
close all
load('F:\真假新闻\script\behav\sub_group.mat');
labels(ID_f)=1;
labels(ID_t)=2;

labels_bin = double(labels == 2);

workpath='.\HMM_results\';
cd(workpath)
dest='.\HMM_results\group_classify';
list=dir('state*');

load('rng_state.mat', 's') % seed 
for i=1:11
    cd([workpath,list(i).name]);
    list2=dir('Sum*.mat');
    load(list2(1).name);
    FO(isnan(FO))=0;
    avg_life(isnan(avg_life))=0;
    features=[FO(1:60,:),avg_life(1:60,:)]; % for listening 
   % features=[FO(61:120,:),avg_life(61:120,:)]; % for recall 
   

        
   % five-fold cross-vidliation
    rng(s); 
    cv = cvpartition(labels_bin, 'KFold', 5);
    accuracy = zeros(cv.NumTestSets, 1);
    n_features=size(features,2);
    weights = zeros(cv.NumTestSets, n_features);
    for k = 1:cv.NumTestSets
        train_idx = training(cv, k);
        test_idx = test(cv, k);     
        X_train = features(train_idx, :);
        y_train = labels_bin(train_idx);
        X_test = features(test_idx, :);
        y_test = labels_bin(test_idx);

        % 使用线性核 SVM
        SVMModel = fitcsvm(X_train, y_train, 'KernelFunction', 'linear', 'Standardize', true); %X_train, y_train
        tmp = abs(SVMModel.Beta);
        %tmp=tmp./sum(tmp);
        weights(k, :)=tmp;
        y_pred = predict(SVMModel, X_test); 
   
        clear tmp
        
                   
        % get accuracy
        accuracy(k) = mean(y_pred == y_test' );
    end
    
    
    result(i)=mean(accuracy);
    mean_weights{i} = mean(weights);

   % sort features
   %[sorted_weights, sorted_idx] = sort(mean_abs_weights, 'descend');
   
    n = 60;
    x = round(max(result(i)) * n);
    p = binocdf(x - 1, n, 0.5);  
    p_values(i) = 1 - p;            

end





