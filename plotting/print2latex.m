function print2latex(filename, figsize, varargin)

set(0, 'defaulttextinterpreter', 'none');

printcmd = sprintf('print(gcf, ''-depsc2'', ''%s'')', filename)
laprint(gcf, filename, 'width', figsize, varargin{:}, 'printcmd', printcmd);