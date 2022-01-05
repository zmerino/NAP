
% class assignment
actual = distributions;
actual.generate_data = false;

actual.dist_name = "Normal";
% file name for actual distribution. "A_" puts at the top of the folder.
actual.filename = sprintf(['A_', char(actual.dist_name),'_Act']);

% creat rndom object
rndom = actual;
rndom.Ns = 2^9;
rndom.randomVSactual = "random";
rndom = dist_list(rndom);
sample = rndom.rndData;

p = [1,0.55,1,0.33,2,ceil(0.0625*rndom.Ns^0.5),40];

filename = "test";
send_file_name = "test2";
min_limit = 0;
max_limit = 1000;


[fail_code,x,SE_pdf,SE_cdf,SE_u,SE_SQR,nBlocks,Blacklist,...
    rndom.Ns,binrndom.Ns, max_LG, sum_LG,T,BRlevel,BR0]...
    = stitch_pdf(sample,filename,send_file_name,min_limit,max_limit,p);

plot(x,SE_pdf)