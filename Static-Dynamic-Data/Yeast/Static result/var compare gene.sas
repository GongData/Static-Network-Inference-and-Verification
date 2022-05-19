
%MACRO compare(file,file_out);
   PROC IMPORT OUT= WORK.data1 
            DATAFILE= "C:\Users\Yinjiao\Desktop\results\&file" 
            DBMS=TAB REPLACE;
        GETNAMES=YES;
        DATAROW=2; 
   RUN;
   Data data2;
        set data1;
        gene_compare=compress(gene1||"--"||gene2||";");
   Run;
   PROC EXPORT DATA= WORK.data2 
            OUTFILE= "C:\Users\Yinjiao\Desktop\&file_out" 
            DBMS=CSV REPLACE;
            PUTNAMES=YES;
   RUN;
%mend compare;
%compare(uniq_t_static_0.3.txt,uniq_t_static_0.3.csv)
%compare(uniq_t_static_0.28.txt,uniq_t_static_0.28.csv)
%compare(uniq_t_static_0.32.txt,uniq_t_static_0.32.csv)
%compare(uniq_t_static_0.34.txt,uniq_t_static_0.34.csv)
