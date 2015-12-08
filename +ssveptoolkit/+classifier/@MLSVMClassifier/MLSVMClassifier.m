classdef MLSVMClassifier < ssveptoolkit.classifier.ClassifierBase
    
    properties (Constant)

    end
    
    properties
       Coding; % Coding design 'onevsall' (default) | 'allpairs' | 'binarycomplete' | 'denserandom' | 'onevsone' | 'ordinal' | 'sparserandom' | 'ternarycomplete' | numeric matrix
       FitPosterior; % Flag indicating whether to transform scores to posterior probabilities false or 0 (default) | true or 1
       Prior % 'empirical' (default) or 'uniform'.  Prior probabilities for each class.
       tmpltSVM;
       models; 
    end
    
    methods (Access = public)
        function MLSVM = MLSVMClassifier(instanceSet, coding, fitPosterior, prior, tmpltsvm)
            %set default parameters
            
            if nargin == 0
                MLSVM.Coding='onevsone';
                MLSVM.FitPosterior='off';
                MLSVM.Prior='empirical';
                MLSVM.tmpltSVM=templateSVM('KernelFunction','linear');
            elseif nargin == 1 
                MLSVM.instanceSet = instanceSet;
                MLSVM.Coding='onevsone';
                MLSVM.FitPosterior='off';
                MLSVM.Prior='empirical';
                MLSVM.tmpltSVM=templateSVM('KernelFunction','linear');              
            elseif nargin == 2 
                MLSVM.instanceSet = instanceSet;
                MLSVM.Coding=coding;
                MLSVM.FitPosterior='off';
                MLSVM.Prior='empirical';
                MLSVM.tmpltSVM=templateSVM('KernelFunction','linear');
            elseif nargin == 3 
                MLSVM.instanceSet = instanceSet;
                MLSVM.Coding = coding;
                MLSVM.FitPosterior=fitPosterior';
                MLSVM.Prior='empirical';
                MLSVM.tmpltSVM=templateSVM('KernelFunction','linear');
            elseif nargin == 4
                MLSVM.instanceSet = instanceSet;
                MLSVM.Coding = coding;
                MLSVM.FitPosterior=fitPosterior;
                MLSVM.Prior=prior;
                MLSVM.tmpltSVM=templateSVM('KernelFunction','linear');
            elseif nargin == 5
                MLSVM.instanceSet = instanceSet;
                MLSVM.Coding = coding;
                MLSVM.FitPosterior=fitPosterior;
                MLSVM.Prior=prior;
                MLSVM.tmpltSVM=tmpltsvm;
            elseif nargin > 5
                MLSVM.instanceSet = instanceSet;
                MLSVM.Coding = coding;
                MLSVM.FitPosterior=fitPosterior;
                MLSVM.Prior=prior;
                MLSVM.tmpltSVM=tmpltsvm;             
            end
        end
        
        function MLSVM = build(MLSVM)
            %clear all from previous calls to "build"
            MLSVM.reset;
            numLabels = MLSVM.instanceSet.getNumLabels;
            uniqueLabels = unique(MLSVM.instanceSet.getLabels);
           
            % ---- Multi-Class ----- %
            instances=MLSVM.instanceSet.instances;
            labels=MLSVM.instanceSet.labels;
            
            %t=templateSVM('KernelFunction','linear');
            MLSVM.models{1}=fitcecoc(instances,labels,'Coding', MLSVM.Coding,'FitPosterior', MLSVM.FitPosterior,'Prior',MLSVM.Prior,'Learners',MLSVM.tmpltSVM); 
            % ---- One (vs) All ----- %
%             for i=1:numLabels
%                 currentLabel = uniqueLabels(i);
%                 labels = zeros(LDA.instanceSet.getNumInstances,1)-1;
%                 labels(LDA.instanceSet.getInstanceIndicesForLabel(currentLabel)) = 1;
%                 instances = sparse(LDA.instanceSet.getInstances);
%                 LDA.models{i} = libsvmtrain(labels,instances, '-t 0 -c 1 -b 1');
%                 if LDA.kernel == LDA.KERNEL_LINEAR;
%                    %store the models in an instance variable
%                    LDA.models{i} = libsvmtrain(labels, instances, sprintf('-t %d -c %f -b 1 -q', LDA.kernel, LDA.cost));
%                 elseif LDA.kernel == LDA.KERNEL_RBF;
%                    LDA.models{i} = libsvmtrain(labels, instances, sprintf('-t %d -c %f -g %f -b 1 -q', LDA.kernel, LDA.cost, LDA.gamma));
%                 else
%                    error('invalid kernel parameter');
%                 end
%             end
        end
        
        function [output, probabilities, ranking] = classifyInstance(MLSVM,instance)
            %input = instance matrix rows = instances, cols = attributes
            %output = predicted class
            %probabilities = probability for predicted class
            %ranking = propabilities for all classes (e.g. to use with mAP)
            
            
            %TODO:should print an error if 'build' has not been called
            numModels = length(MLSVM.models);
            [numinstance, ~] = size(instance);
            %scores = zeros(numModels,numinstance);
             
            % ---- Multi-class ----- %
            [label,scores,loss] = predict(MLSVM.models{1},instance); 
            
            % ---- One (vs) All -----%
%              for i=1:numModels
%                  %predict using the stored models
%                  [label,score,cost] = predict(LDA.models{i},instance);
%                  %libsvmpredict(eye(numinstance,1),instance, LSVM.models{i},'-b 1 -q');
%                 %store probability for each class
%                 scores(i,:) = score(:,1);
%             end

             output = zeros(numinstance,1);
             probabilities = zeros(numinstance,1);
             %we need these for ranking metrics (e.g. mAP)
             ranking = scores;
             for i=1:numinstance
                 %select the class with the highest probability
                 [prob, idx] = max(scores(i,:));
                 uniqueLabels = unique(MLSVM.instanceSet.getLabels);
                 %output the label with highest probability
                 output(i,1) = uniqueLabels(idx);
                 %return the probability for the output label
                 probabilities(i,1) = prob;
             end
        end
        
        function MLSVM = reset(MLSVM)
            %delete all stored models
            MLSVM.models = {};
        end
        
        function configInfo = getConfigInfo(MLSVM)
            %configInfo=sprintf('MLSVMClassifier\tCoding:%s\tFitPosterior:%s\tPrior:%s\templateSVM-Kernel:%s\n', MLSVM.Coding, MLSVM.FitPosterior,MLSVM.Prior,MLSVM.tmpltSVM);
            configInfo=sprintf('MLSVMClassifier\tCoding:%s\tFitPosterior:%s\tPrior:%s\n', MLSVM.Coding, MLSVM.FitPosterior,MLSVM.Prior);
            disp(MLSVM.tmpltSVM)
            %configInfo = 'MLSVMClassifier (Config info not supported yet)';
        end
                
    end
end

