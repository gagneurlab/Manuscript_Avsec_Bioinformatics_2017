

## TODO - add protein-labels

#'
#' ## GAM vs no position

ggplot(dt[task != "branch-point"], aes(x = no_pos, y = gam,
                                  color = task, shape = signif.gam__no_pos, alpha = signif.gam__no_pos)) +
  geom_point() +
  geom_abline() +
  scale_shape_manual(values=c(19,1)) +
  scale_alpha_manual(values=c(0.6, 1)) + 
  facet_wrap(~metric+task, scales="free", ncol=2) +
  legend_top

## on how many proteins do we see an improvement
dt[task != "branch-point", .(signif_better = sum(gam > no_pos& signif.gam__no_pos),
       signif_worse = sum(gam < no_pos & signif.gam__no_pos),
       ns = sum(!signif.gam__no_pos)
    ) , by = .(metric, task)]

## o and x ?

##' gam vs relu
ggplot(dt[order(task, decreasing = TRUE)], aes(x = relu, y = gam,
               color = task, alpha = signif.gam__relu, shape = signif.gam__relu)) +
  geom_point() +
  geom_abline() +
  facet_wrap(~metric, scales="free") +
  scale_shape_manual(values=c(19,1)) +
  scale_alpha_manual(values=c(0.6, 1)) + 
  legend_top

ggplot(dt[order(task, decreasing = TRUE)], aes(x = relu, y = gam,
               color = signif.gam__relu, shape = signif.gam__relu)) +
  geom_point() +
  geom_abline() +
  facet_wrap(~metric + task, scales="free") +
  scale_shape_manual(values=c(19,1)) +
  scale_alpha_manual(values=c(0.6, 1)) + 
  legend_top


dt[, .(signif_better = sum(gam > relu & signif.gam__relu),
       signif_worse = sum(gam < relu & signif.gam__relu),
       ns = sum(!signif.gam__relu)
    ) , by = .(metric, task)]
