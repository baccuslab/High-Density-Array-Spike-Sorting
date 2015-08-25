function [ok] = opt_cbar(cax)

if ( nargin ~= 1 )
  
  cax = caxis;
  
end
oax = gca;
p=get(gca,'position');
p(1)=p(1)+p(3)+.01;
p(3)=.02;
hll=axes('position',p);
axes(hll)
scale=[1:length(colormap)]';
xlim=[0 1];
ylim=[cax(1) cax(2)*1.00001];
image(xlim,ylim,scale)
set(hll,'ydir','normal','ylim',ylim,'box','on','xtick',[],'yaxislocation','right');

set(hll,'ticklength',[0 0],'yscale','linear')
axis xy
axes(oax)
