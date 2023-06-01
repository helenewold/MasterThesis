# Capon difference plots
When both the USTB minimum_variance_capon and the getCapon postprocess has
been implemented, the difference between the results may be shown by:
$$
diff\_result = abs\left( \frac{A}{max(A(:))} - \frac{B}{max(B(:))} \right)
$$
By applying this with the result from the minimum_variance_capon as A and 
getCapon as B, we get the following result

<p float="left">
  <img src=../Figures/Capon_Section/P4_FI_121444_45mm_focus/diff_P4_FI_121444_45mm_focus.png width=40% />
</p>

Here, it is possible to see a zebra-looking pattern in the depth direction.
The most likely explanation to this is due to numerical differences in the 
process of choosing elements and subarrays in the subarray averaging method.
The two methods may choose the new M, i.e. the number of active elements,
slightly differently, leading to differences in the numerical results.

When visually analyzing the two plots, the differences are visually
insignificant. If one were to look for a very long time, a slight difference
could be possible to notice.

This zebra pattern shows us that this method of comparison is not the most
robust visualization of differences between the two methods, as the
preferred result is an entirely black image, assuming the two methods are
equally good. A slight numerical difference is expected to an extent, as
Matlab may encounter different rounding methods due to small deviations in
the two calculations.

Either way, throughout the thesis, this kind of result
comparison is not useful, as the methods that are implemented are not the 
same, and the results are exptected to be different. Other methods of
comparison, such as resolution or contrast comparisons, are more relevant.
