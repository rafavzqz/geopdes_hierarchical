function bnd_indices = get_boundary_indices (iside, size_dir, indices)

  ndim = numel (size_dir);
  ind2 = ceil (iside/2);
  ind = setdiff (1:ndim, ind2);

  if (mod(iside,2) == 1)
    boundary_ind = 1;
  else
    boundary_ind = size_dir (ind2);
  end

  indsub = cell (1, ndim);
  [indsub{:}] = ind2sub (size_dir, indices);
  aux = find (indsub{ind2} == boundary_ind);
  ppp = cellfun (@(x) x(aux), indsub(ind), 'UniformOutput', false);
  bnd_indices = sub2ind ([size_dir(ind), 1], ppp{:});

end