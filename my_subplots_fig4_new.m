function my_subplots_fig4_new(nA,S,v,w,TCcorr,SMcorr,rTC,rSM,K)
    N = size(rTC{1},1);
    axis off
    set(gca,'Units','normalized','Position',[0 0 1 1]);
    vec = 0.0175:1/nA:1; 
    text(0.10,0.02, '(a)','Color','b','FontSize',14)
    text(0.26,0.02, '(b)','Color','b','FontSize',14)
    text(0.42,0.02, '(c)','Color','b','FontSize',14)
    text(0.59,0.02, '(d)','Color','b','FontSize',14)
    text(0.76,0.02, '(e)','Color','b','FontSize',14)
    text(0.93,0.02, '(f)','Color','b','FontSize',14)

    text(0.16,0.05, ['$\bar{\rho}$ =' num2str(round(sum(SMcorr(2,:))/K,2),'%0.2f') ],'Color','k','FontSize',14, 'Interpreter','latex')
    text(0.33,0.05, ['$\bar{\rho}$ =' num2str(round(sum(SMcorr(3,:))/K,2),'%0.2f') ],'Color','k','FontSize',14, 'Interpreter','latex')
    text(0.50,0.05, ['$\bar{\rho}$ =' num2str(round(sum(SMcorr(4,:))/K,2),'%0.2f') ],'Color','k','FontSize',14, 'Interpreter','latex')
    text(0.67,0.05, ['$\bar{\rho}$ =' num2str(round(sum(SMcorr(5,:))/K,2),'%0.2f') ],'Color','k','FontSize',14, 'Interpreter','latex')
    text(0.84,0.05, ['$\bar{\rho}$ =' num2str(round(sum(SMcorr(6,:))/K,2),'%0.2f') ],'Color','k','FontSize',14, 'Interpreter','latex')

    text(0.25,0.05, ['$\bar{\gamma}$ =' num2str(round(sum(TCcorr(2,:))/K,2),'%0.2f') ],'Color','k','FontSize',14, 'Interpreter','latex')
    text(0.41,0.05, ['$\bar{\gamma}$ =' num2str(round(sum(TCcorr(3,:))/K,2),'%0.2f') ],'Color','k','FontSize',14, 'Interpreter','latex')
    text(0.58,0.05, ['$\bar{\gamma}$ =' num2str(round(sum(TCcorr(4,:))/K,2),'%0.2f') ],'Color','k','FontSize',14, 'Interpreter','latex')
    text(0.75,0.05, ['$\bar{\gamma}$ =' num2str(round(sum(TCcorr(5,:))/K,2),'%0.2f') ],'Color','k','FontSize',14, 'Interpreter','latex')
    text(0.92,0.05, ['$\bar{\gamma}$ =' num2str(round(sum(TCcorr(6,:))/K,2),'%0.2f') ],'Color','k','FontSize',14, 'Interpreter','latex')
    

    for i =1:nA 
        for j=1:S
            ihs = 0.00;  %initial_horizontal_shift
            ivs = 1.00;  %initial_vertical_shift
            shz = 0.05; %subplot_horizontal_size
            svs = S;  %subplot_vertical_size (more the better)
            hs  = 0.725;  %horizontal shift of subplots
            nR  = S+0.5; %No. of rows
            shifter = 0.23;

            if i==1
            %%
            hax=axes();
            imagesc(flipdim(reshape(abs(zscore(rSM{i}(j,:))),v,w),1)); 
            newPos=[hs*(mod(j-1,1)+ihs+0.00),   (1-1/nR)-(1/nR)*(fix((j-1)/1)+ivs-1),   shz,   1/svs];
            set(gca,'outer',newPos), 

            hax=axes();
            plot(zscore(rTC{i}(:,j)));  axis([0 N -3 3]);
            newPos=[hs*(mod(j-1,1)+ihs+0.04),   (1-1/nR)-(1/nR)*(fix((j-1)/1)+ivs-1),   3*shz,   1/svs];
            set(gca,'outer',newPos), 

            else

            %%
            zscore_rxSM = abs(zscore(rSM{i}(j,:)));
            hax=axes();
            imagesc(flipdim(reshape(zscore_rxSM,v,w),1));  colormap('hot')
            newPos=[hs*(mod(j-1,1)+ihs+shifter*(i-1)),   (1-1/nR)-(1/nR)*(fix((j-1)/1)+ivs-1),   shz,   1/svs];
            set(gca,'outer',newPos), 
            set(gca,'XTickLabel','')
            set(gca,'YTickLabel','')
            xlabel(['\rho',' = ',num2str(round(SMcorr(i,j),2))],'color','r','FontSize',12)

            hax=axes();
            plot(zscore(rTC{i}(:,j))); axis([0 N -3 3]);
            newPos=[hs*(mod(j-1,1)+ihs+shifter*(i-1)+0.035),   (1-1/nR)-(1/nR)*(fix((j-1)/1)+ivs-1),   3*shz,   1/svs];
            set(gca,'outer',newPos),
            set(gca,'XTickLabel','')
            set(gca,'YTickLabel','')
            xlabel(['\gamma',' = ',num2str(round(TCcorr(i,j),2))],'color','r','FontSize',12)
            end
        end
    end
