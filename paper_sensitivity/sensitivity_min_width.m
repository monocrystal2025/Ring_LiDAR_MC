function [minimum,gridMinimum] = sensitivity_min_width(width,probability,target)
%SENSITIVITY_MIN_WIDTH First raw crossing, interpolated only between neighbors.
% Columns are separate response curves. No monotone envelope or smoothing is
% applied: a larger divergence is not assumed to preserve detection success.
arguments
    width (:,1) double {mustBeFinite,mustBePositive}
    probability double {mustBeFinite,mustBeNonnegative}
    target (1,1) double {mustBePositive}
end
assert(issorted(width,'strictascend') && numel(width)>=2,'Sensitivity:Grid','Widths must strictly increase.');
assert(size(probability,1)==numel(width) && all(probability<=1,'all') && target<=1, ...
    'Sensitivity:Probability','Invalid probability curve.');
[attained,index]=max(probability>=target,[],1);
previous=max(1,index-1);
columns=1:size(probability,2);
lo=probability(sub2ind(size(probability),previous,columns));
hi=probability(sub2ind(size(probability),index,columns));
minimum=width(previous).'+(width(index)-width(previous)).'.*(target-lo)./max(hi-lo,eps);
minimum(index==1)=width(1);
minimum(~attained)=nan;
gridMinimum=width(index).';gridMinimum(~attained)=nan;
end
