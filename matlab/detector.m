function [b_hat, nerr] = detector(y, b)
% DETECTOR  Threshold detector and error counter (Project 12.7, block 5).
%   [b_hat, nerr] = detector(y, b) decides 1 when y > 0 and 0 otherwise,
%   then counts the decisions that disagree with the transmitted bits b.
%   (y == 0 has probability zero for continuous noise.)
    b_hat = double(y > 0);
    nerr  = sum(b_hat ~= b);
end
