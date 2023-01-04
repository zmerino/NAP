
% sxb2 = load("sample_block_2.txt");
% tempStruc.lowBound = min(sxb2);
% tempStruc.highBound = max(sxb2);
% [~, x2, pdf2, cdf2, ~,~] = EstimatePDF(sxb2,tempStruc);
% 
% sxb3 = load("sample_block_3.txt");
% tempStruc.lowBound = min(sxb3);
% tempStruc.highBound = max(sxb3);
% [~, x3, pdf3, cdf3, ~,~] = EstimatePDF(sxb3,tempStruc);


sxb2 = load("sample_block_11.txt");
tempStruc.lowBound = min(sxb2);
tempStruc.highBound = max(sxb2);
[~, x2, pdf2, cdf2, ~,~] = EstimatePDF(sxb2,tempStruc);

sxb3 = load("sample_block_12.txt");
tempStruc.lowBound = min(sxb3);
tempStruc.highBound = max(sxb3);
[~, x3, pdf3, cdf3, ~,~] = EstimatePDF(sxb3,tempStruc);

% sxb2 = load("sample_block_2.txt");
% [~, x2, pdf2, cdf2, ~,~] = EstimatePDF(sxb2);
% 
% sxb3 = load("sample_block_3.txt");
% [~, x3, pdf3, cdf3, ~,~] = EstimatePDF(sxb3);
% 
% sxb2 = load("sample_block_11.txt");
% [~, x2, pdf2, cdf2, ~,~] = EstimatePDF(sxb2);
% 
% sxb3 = load("sample_block_12.txt");
% [~, x3, pdf3, cdf3, ~,~] = EstimatePDF(sxb3);

disp(['Input Block 2: min ', num2str(min(sxb2)), ' max ', num2str(max(sxb2))])
disp(['Input Block 2: min ', num2str(min(sxb3)), ' max ', num2str(max(sxb3))])
disp(['Output Block 2: min ', num2str(min(x2)), ' max ', num2str(max(x2))])
disp(['Output Block 2: min ', num2str(min(x3)), ' max ', num2str(max(x3))])

figure()
hold on;
plot([min(x2),max(x2)],[2,2],'-ro')
plot([min(sxb2),max(sxb2)],[2.5,2.5],'--rs')
plot([min(x3),max(x3)],[3,3], '-bo')
plot([min(sxb3),max(sxb3)],[3.5,3.5], '--bs')






