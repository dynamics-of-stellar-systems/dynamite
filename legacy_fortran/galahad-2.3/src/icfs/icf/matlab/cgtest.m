clear

prob_list  = char('1138bus','bcsstk08', 'bcsstk09', 'bcsstk10', 'bcsstk11', ...
    'bcsstk18', 'bcsstk19') ;

for i = 1:size(prob_list,1) 
   load(['../tprobs/', deblank(prob_list(i,:))]) ;
   
   n = size(A,1) ;  
   Adiag = full(diag(A)) ;   lA = tril(A) - spdiags(Adiag, 0, n, n) ;
   maxiter = max(n, 100) ;
   
   p = 5 ;
   [lL, Ldiag] = icfmex(n, lA, Adiag, p) ;
   
   rtol = 1.0e-6 ;
   b = A*ones(n,1) ;

   [x, cgiter, info] = cgmex(n, lA, Adiag, lL, Ldiag, b, maxiter, rtol) ;
   fprintf('nnz in L = %d\n', nnz(lL)) ;
   fprintf('cgiter = %d info =  %d relative error %.12f\n', cgiter, info, ...
	   norm(A*x-b)/norm(b)) ;
end
