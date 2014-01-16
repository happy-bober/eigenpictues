#The procedure eigenpicture(A,n,shape) generates an eigenpicture of the 2x2 matrix A using n unit circle vectors and plots the shape created by applying A to the vectors if the value of the shape parameter is true. On the picture the eigenvectors of A and their transformations are displayed in black (if they are real)  and the eigenvectors of 
#A^T.A;
# are displayed in blue (if they are real).Note that n and shape are optional arguments, their default values are n=40 and shape=false.

eigenpicture:=proc(A,n::integer :=40,shape::truefalse := false)  
#A is the matrix generating the eigenpicture; n is the number of vectors plotted, 40 by default; shape=true plots an outline of the eigenpicture
  local m, e_value,e_vectors,ATA,s_value,s_vectors,base_vectors,trans_vectors,base_plot,trans_plot,e_plot1,e_plot2,s_plot1,s_plot2,poly:
  m:=2:                                                        #assumed dimension of the matrix (2)
  e_value,e_vectors:=LinearAlgebra[Eigenvectors](A);           #get eigenvectors of the matrix
  if e_value[1]=e_value[2] then:
    e_vectors:=[seq(e_vectors[..,i], i=1..m)];
  else:
    e_vectors:=[seq(e_vectors[..,i]*(1/LinearAlgebra[VectorNorm](e_vectors[..,i],Euclidean)), i=1..m)]; 
                                                               #generate eigenvectors and normalise them unless the eigenvalues are both 0, then the second eigenvector is of length 0 and thus we can't normalise it
  end if:
  ATA:=LinearAlgebra[Transpose](A).A;

  s_value,s_vectors:=LinearAlgebra[Eigenvectors](ATA);         #get eigenvectors of A^T.A
  if s_value[1]=s_value[2] then:
    s_vectors:=[seq(s_vectors[..,i], i=1..m)];
  else:
    s_vectors:=[seq(s_vectors[..,i]*(1/LinearAlgebra[VectorNorm](s_vectors[..,i],Euclidean)), i=1..m)]; 
                                                               #generate eigenvectors of A^T.A and normalise them unless the eigenvalues are both 0, then the second eigenvector is of length 0 and thus we can't normalise it
  end if:
  base_vectors:=[seq(<cos(2*k*Pi/n),sin(2*k*Pi/n)>, k=1..n)]:  #generate unit circle vectors
  trans_vectors:=map(x->A.x,base_vectors):                     #apply A to the unit circle vectors
  base_plot:=[seq(plottools[arrow]([0,0], base_vectors[k],0.1,0.1,0.01,arrow,color=ColorTools:-Color("HSV",[k/n,0.5,1])),k=1..n)]:
                                                               #generate unit circle arrow plots
  trans_plot:=[seq(plottools[arrow](base_vectors[k],trans_vectors[k],0.1,0.1,0.01,arrow,color=ColorTools:-Color("HSV",[k/n,1,1])),k=1..n)]:
                                                               #generate transformed unit circle arrow plots
  e_plot1:=[seq(plottools[arrow]([0,0],e_vectors[k],0.5,0.1,0.1,arrow),k=1..m)]:                                                                                            #generate eigenvector arrow plots
  e_plot2:=[seq(plottools[arrow](e_vectors[k],A.e_vectors[k],0.5,0.1,0.01,arrow),k=1..m)]:                           s_plot1:=[seq(plottools[arrow]([0,0],s_vectors[k],0.3,0.1,0.1,arrow,color="Blue"),k=1..m)]:
                                                               #generate A^T.A eigenvector arrow plots
  s_plot2:=[seq(plottools[arrow](s_vectors[k],A.s_vectors[k],0.3,0.1,0.01,arrow,color="Blue"),k=1..m)]:
  if shape then                                                #create the outline plot
    poly:=plots[polygonplot](map(x->x[1],[seq(base_vectors[k]+trans_vectors[k],k=1..n)]),map(y->y[2],[seq(base_vectors[k]+trans_vectors[k],k=1..n)]),  style=line):
    plots[display]([poly,op(base_plot),op(trans_plot),op(e_plot1),op(e_plot2),op(s_plot1),op(s_plot2)],scaling=constrained);
  else:                                                        #plot everything on one picture
    plots[display](base_plot,trans_plot,e_plot1,e_plot2,s_plot1,s_plot2,scaling=constrained);
  end if:
end proc:

#Here is a script to animate an eigenpicture of cos(t)*[2  1/2]
#                                                      [-1  2 ] 
plots[animate](eigenpicture,[Matrix([[2*cos(t),1/2*cos(t)],[-cos(t),2*cos(t)]]),80,true], t=0..2*Pi,scaling=constrained,frames=100);
