%Multi-view Clustering using a flexible and optimal multi-graph fusion method
%School of Computer Science and Information, AnHui Polytechnic University, Wuhu, 241000, Anhui, China
%X:data      
%gt:number of clusters
%alf,gamma and beita are parameters
%nCluster represents the number of classes
%s represents a custom column sum
%
function [P,Zx ]=Zv_and_Zx(X,gt,alf,beita,gamma,nCluster,s )
%∑||Xv-XvZv||^2+alf||Zv||^2+beita*||Zv-Z*||^2+gammaTr(P'LP)
%  s.t. 1'Z*=s1'   P'P=I
    iter=0;
    iterMax=35;
    m=size(X,1);
    n=size(X{1},2);
    c=length(unique(gt));
    Convalue=[];
    P=zeros(n,c);
    Zx=zeros(n);
    for i=1:m
       Z{i}=zeros(n);
       K{i}=X{i}'*X{i};
    end
    D=zeros(n);
    while iter<iterMax
        value=0;
     % fprintf('   iter=%d',iter);
      %更新Zv
      iter=iter+1;
        for i=1:m
           % Z{i}=inv(2*K{i}+2*(alf+beita)*eye(n))*(2*K{i}+2*beita*Zx) ;  
            Z{i}=(2*K{i}+2*(alf+beita)*eye(n))\(2*K{i}+2*beita*Zx) ;
            
        end
      %更新Z*
      temp=zeros(n);
      for i=1:m
         temp=temp+Z{i}; 
      end
          Zx=(2*beita*temp+gamma*(P*P'))/(2*m*beita);      
          Zx=SimplexProj( Zx'/s);  
          Zx= Zx'*s;
         
%        Zx= solver_BCLS_closedForm(Zx/s);
%        Zx=Zx*s;     
      %更新P
     
      S=(Zx+Zx')/2;
      D=diag(sum(S));
      L=D-S;
%       [uN,sN,vN] = svd(L);
%       P=uN(:,n-c+1:n);
   [P,~,~]=eig1(L,c,0);
     
   
    
  % [nmi(iter),acc(iter),f(iter),purity(iter)]=zhixing_Kmeans(P,nCluster,gt);
%     for i=1:m
%        tt=K{i}-2*K{i}*Z{i}+Z{i}*K{i}*Z{i};
%        ttt=Z{i}-Zx;
%        value=value+trace(tt)+alf*trace(Z{i}*Z{i}')+beita*trace(ttt'*ttt);   
%     end
%     value=value+gamma*trace(P'*L*P);
%     Convalue(iter)=value;
    end
    
 %result=[nmi;acc;f;purity];
 
end