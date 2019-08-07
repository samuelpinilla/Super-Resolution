function myfiles(name)

  f=matlab.codetools.requiredFilesAndProducts(name);

  for n=1:length(f)
    fprintf('%s\n',f{n});
  end   

end