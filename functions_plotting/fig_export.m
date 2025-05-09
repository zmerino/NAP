function fig_export(handle,fig_name,fig_style)


style = hgexport('readstyle',fig_style);
style.Format = 'png';
hgexport(handle,fig_name,style);
style.Format = 'fig';
hgexport(handle,fig_name,style);
end
