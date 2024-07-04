function [zero_rows,zero_cols] = get_zero_rows_and_columns(A,method)
    % Gets the rows and columns of matrix A, designed for A to be sparse
    % (but may also work well for full matrices, who knows)

    switch method
        case 'all(A == 0)' % This method has been shown to be slow
            zero_cols = all(A == 0,1); % Sparse logicals are probably a bad way to store these
            zero_rows = all(A == 0,2);
        case 'find(A) - traditional indices'
            [nz_rows_all,nz_cols_all] = find(A);
            zero_rows = setdiff(1:size(A,1),nz_rows_all);
            zero_cols = setdiff(1:size(A,2),nz_cols_all);
            % nz_rows = unique(nz_rows_all);
            % nz_cols = unique(nz_cols_all);
            % zero_rows = setdiff(1:size(Kr,1),nz_rows);
            % zero_cols = setdiff(1:size(Kr,1),nz_cols);
        case 'find(A) - logical indices'
            [nz_rows_all,nz_cols_all] = find(A);
            nz_rows = unique(nz_rows_all);
            nz_cols = unique(nz_cols_all);
            zero_rows = true(size(A,1),1);
            zero_cols = true(size(A,2),1);
            
            zero_rows(nz_rows) = false;
            zero_cols(nz_cols) = false;
        otherwise
            error('method not recognized')
    end
end
