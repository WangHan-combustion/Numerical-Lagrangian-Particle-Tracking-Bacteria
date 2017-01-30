function [TD]=TaylorDispersion(lcx,lcy,mv)
tic
NT=size(lcx,1);
Tau=[1:2:8 fix(10.^[1:0.1:log10(NT-100)])];
NTau=length(Tau);
BINA=[-6:1e-2:3];
NB=length(BINA);
JPDFraw=zeros(NTau,NB);
JPDFrem=zeros(NTau,NB);
Q=0:0.25:8;
NQ=length(Q);
SPDraw=zeros(NTau,NQ);
SPDrem=zeros(NTau,NQ);
for i=1:length(Tau)
    tmpx=lcx(Tau(i)+1:end,:)-lcx(1:end-Tau(i),:);
    tmpy=lcy(Tau(i)+1:end,:)-lcy(1:end-Tau(i),:);
    dr=tmpx.^2+tmpy.^2;
    dr=reshape(dr,[],1);
    dr=sqrt(dr);
    dr=dr(dr>0);
    tmp1=dr.^Q(1);
    tmp2=dr.^0.25;
    [pd,~]=fastpdf(log10(dr),BINA,'not');
    JPDFraw(i,:)=pd;
    for j=1:NQ
        SPDraw(i,j)=mean(tmp1);
        tmp1=tmp1.*tmp2;
    end
    
    tmpx=lcx(Tau(i)+1:end,:)-lcx(1:end-Tau(i),:)-mv(1)*Tau(i)/40;
    tmpy=lcy(Tau(i)+1:end,:)-lcy(1:end-Tau(i),:)-mv(2)*Tau(i)/40;
    dr=tmpx.^2+tmpy.^2;
    dr=reshape(dr,[],1);
    dr=sqrt(dr);
    dr=dr(dr>0);
    tmp1=dr.^Q(1);
    tmp2=dr.^0.25;
    [pd,~]=fastpdf(log10(dr),BINA,'not');
    JPDFrem(i,:)=pd;
    for j=1:NQ
        SPDrem(i,j)=mean(tmp1);
        tmp1=tmp1.*tmp2;
    end
end
TD.SPDraw=SPDraw;
TD.SPDrem=SPDrem;

TD.Tau=Tau;
TD.Q=Q;
TD.JPDFraw=JPDFraw;
TD.JPDFrem=JPDFrem;
TD.BINA=BINA;