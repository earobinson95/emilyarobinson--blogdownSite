x = seq(1, 10000)
y = 1 - (1/(x^2))

value <- NA
for (i in 1:10000){
value[i] <- sum(y[1:i])
}
value
plot(x, value)
