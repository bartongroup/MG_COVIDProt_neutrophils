targets_profiles <- function() {
  profiles <- list(
    tar_target(time_prof, make_time_profiles(covid, treat = "Brensocatib")),
    tar_target(time_prof_cl6, cluster_time_profiles(time_prof, n_clust = 6)),
    tar_target(fig_profiles_cl6, plot_cluster_profiles(time_prof_cl6)),
    tar_target(fig_prof_heatmap, plot_prof_heatmap(time_prof)),
    tar_target(fig_gradients_eigen, plot_gradients(time_prof))
  )
  
  c(
    profiles
  )
}