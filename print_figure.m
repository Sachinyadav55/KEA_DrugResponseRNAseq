function print_figure(h,papersize, output_file,opts)
    set(h, 'PaperPosition', [0 0 papersize]); 
    set(h, 'PaperSize', papersize); 
    
    print(h, opts, output_file);
end