function [X,Y,ncon, n_nodes, n_elements] = GenerateCustomToolMesh(Xb,Yb, L_ref, nx_ref, nx_coarse, ny)
        
%This represents the difention sof the top edge. We split the top edge
%depending on the user defined "refined zone"

        %X points
        xL = Xb(4);%Top leftest point on the tool geometry
        xR = Xb(3);%Top rightest point on the tool geometry
        
        %Y points
        yL = Yb(4);%Top leftest point on the tool geometry
        yR = Yb(3);%Top rightest point on the tool geometry

        L_total = xR - xL;

        %split point
        x_split = xL + L_ref;
        
        %refined region
        x_top1 = linspace(xL, x_split, nx_ref);
        
        %coarse region
        x_top2 = linspace(x_split, xR, nx_coarse);
        
        %combine them and remove duplicate
        x_top = [x_top1, x_top2(2:end)];
        
        % interpolate y to make a straight line
        y_top = yL + (yR - yL) * (x_top - xL) / L_total;
        
        nx = length(x_top);
        
        %This is defining the bottum edge, same method as before
        xL_bot = Xb(1);
        xR_bot = Xb(2);
        
        yL_bot = Yb(1);
        yR_bot = Yb(2);
        
        x_bot = xL_bot + (xR_bot - xL_bot) * (x_top - xL) / L_total;
        y_bot = yL_bot + (yR_bot - yL_bot) * (x_top - xL) / L_total;

        %This is the section where we interpolate the mesh
        t_vals = linspace(0,1,ny);
        
        X = [];
        Y = [];
        %Transfinite interpolation in the thickness direction of the tool
        
        for j = 1:ny
            %Interior node rows are computed as a weighted blend
            %Nodes are distributed smoothly between the top and bottom
            %edges with veritcal columns corresponding to top and bottom
            %points. Boundary node positions are preserved while the
            %interior nodes are interpolated alebraically
            tau = t_vals(j);
        
            x_row = (1-tau)*x_top + tau*x_bot;
            y_row = (1-tau)*y_top + tau*y_bot;
        
            X = [X; x_row'];
            Y = [Y; y_row'];
        end
        
        %Generating the connectivity bewtween the nodes witha double loop
        ncon = [];
        %Explantion of the loop can be found in the references section of
        %the report. But the consensus is that for each quadrilateral cell
        %formed by four neighbouring nodes, two triangles are created by
        %splitting along the diagonal
        for j = 1:(ny-1)
            for i = 1:(nx-1)
                n1 = (j-1)*nx + i;
                n2 = n1 + 1;
                n3 = n1 + nx;
                n4 = n3 + 1;
        
                ncon = [ncon; n1 n2 n3];
                ncon = [ncon; n2 n4 n3];
            end
        end
        
        n_nodes = length(X);
        n_elements = size(ncon,1);




end