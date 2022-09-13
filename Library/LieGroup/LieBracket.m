function adxy = LieBracket(x,y)
% Lie Bracket Operator on Algebra of a Matrix Lie Group
% From Lie Algebras (twists) to Lie Algebra (twists)
% (N x N) g x g -> (N x N) g
    
    adxy = x*y - y*x;
end