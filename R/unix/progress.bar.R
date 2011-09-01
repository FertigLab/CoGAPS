updatePB <- function(start.iter, end.iter, adapting)
{
  tcltk:::tkProgressBar(title=ifelse(adapting, "Adapting","Updating"),
                        label="Iteration 0", min = start.iter, max=end.iter,
                        initial=start.iter)
}

setPB <- function(pb, iter)
{
  tcltk:::setTkProgressBar(pb, iter, label=paste("Iteration",iter))
}
