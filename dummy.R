library("ggplot2")

# Create a dummy plot
p <- ggplot(mtcars, aes(x=wt, y=mpg)) + geom_point()

# SHow the plot
print(p)

