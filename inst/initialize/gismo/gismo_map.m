% GISMO_MAP: internal function to obtain the mapping of a G+smo geometry,
% and its first and second derivatives. 
%
% function F = gismo_map (pts, in, 1)
%          jac = gismo_map (pts, in, 2)
%          [F, jac] = gismo_map (pts, in, 2)
%          hess = gismo_map (pts, in, 3)
%          [F, hess] = gismo_map (pts, in, 2)
%          [F, jac, hess] = gismo_map (pts, in, 2)
%
% INPUT:
%   pts: cell or matrix, points at which the mapping is evaluated
%   in: gsTHBSpline or gsTensorBspline object coming from G+smo
%   der: 0 to get the image of pts by the geometric mapping of in
%        1 go get the first derivative of pts
%        2 to get the second derivative
%
% OUTPUT: 
%   F: [rdim x npts] G+smo geometry evaluated at pts
%   jac: [rdim x ndim x npts] jacobian of the G+smo geometry evaluated at pts
%   hess: [rdim x ndim x ndim x npts] hessian of the G+smo geometry evaluated at pts
%
% Copyright (C) 2018 Ondine Chanon


function varargout = gismo_map (pts, in, der)
  if ( iscell (pts) ) 
      ndim = length (pts);
      npts = prod (cellfun (@length,pts));
      pts = cartesian_product_from_cell (pts);
  else
      [ndim, npts] = size (pts);
  end
  
  if (der == 0 || nargout > 1)
    F = in.eval(pts); % dimension rdim x npts
    varargout{1} = F;
  end
  if (der == 1 || nargout > 2)
    % g+smo dim: rdim x (ndim x npts)
    % geopdes dim: rdim x ndim x npts
    jac = reshape (in.jacobian(pts), [], ndim, npts); 
    if nargout == 1
        varargout{1} = jac;
    else
        varargout{2} = jac;
    end
  end
  if (der == 2)
    rdim = in.geoDim;
    % g+smo dim: rdim, (ndim x ndim) x npts
    % geopdes dim: rdim x ndim x ndim x npts
    hess = zeros (rdim, ndim, ndim, npts);
    if ndim ~= 1
      for dir = 1:rdim
        hess(dir,:,:,:) = reshape (in.hess(pts, dir), ndim, ndim, npts);
      end
    else
      % TODO Hessian in G+smo
      warning('Hessian not implemented yet in G+smo for geometries with parametric dimension 1')
    end
    varargout{nargout} = hess;
  end
end


function pts_aux = cartesian_product_from_cell (pts)
  % create cartesian product points from cell information
  s = cellfun (@length, pts);
  s_cell = cell (length(s), 1);
  for ii = 1:length(s)
      s_cell{ii} = 1:s(ii);
  end
  x = cell (1, numel (s_cell));
  [x{:}] = ndgrid (s_cell{:});
  pts_aux = [];
  for ii = 1:length(s)
      pts_aux = [pts_aux ; reshape( pts{ii}(x{ii}), 1, [] )];
  end
end