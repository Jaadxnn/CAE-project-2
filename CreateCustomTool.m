function [Xb, Yb] = CreateCustomTool (Wt, Wb, H)
            % Check validity
   if Wt < Wb
      error('Top width must be greater than the bottom width');
   end
        
    dx = (Wb - Wt) / 2;

    % Define vertices (bottom on y = 0, centered at x = 0)
    X = [-Wb/2; Wb/2; Wb/2 - dx; -Wb/2 + dx; -Wb/2];
    Y = [0; 0; H; H; 0];
    Xb = X + Wt/2;
    Yb = Y - H;
end

