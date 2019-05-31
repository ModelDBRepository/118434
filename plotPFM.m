
function plotPFM(start,rows,columns,X,Y,HPC,threshold)

	figure('Name','Place Field Maps','NumberTitle','off')
	k=0;
	for i=1:rows
		for j=1:columns
			k=k+1;
			subplot(rows,columns,k);
            caxis([threshold 1]);
            axis off
            hold on
            mesh(X,Y,HPC(:,:,start)), view([0 90]), shading interp
			start=start+1;
		end;
	end;