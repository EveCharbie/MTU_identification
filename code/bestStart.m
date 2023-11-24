function [x0] = bestStart(a,lUMT,BestStartMatrix)
% a : muscle activation (double)
% lUMT : UMT length (double)
% BestStartMatrix : maxtrix of the value of tendonforce(N), muscleforce(N),
% tendonlengthning(m), fibrelength(m), pennationangle(rad)
vectSize = size(BestStartMatrix,1) ; 
vect = ones(vectSize,1) ;

err = sum(abs(BestStartMatrix(:,1:2) - [vect.*a ,vect.*lUMT]),2) ; % error vect 
[~,idx] = min(err) ; % index of minimal value 

x0 = BestStartMatrix(idx,3:end); 
end
