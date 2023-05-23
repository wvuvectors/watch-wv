df_campus <- df_watch %>% filter(type == "higher ed dorm") %>% 
	group_by(day) %>% 
	arrange(day) %>% 
	summarize(day_sum = sum(n1n2), day_trend = sum(n1n2.day5.mean))

p1 <- ggplot(df_campus) + labs(y = "", x = "") + 
	scale_y_continuous(labels = comma) + 
	scale_x_date(breaks = "1 month", labels = format_dates) + 
	scale_color_manual(name = "Target", values = TARGET_COLORS, labels = c("n1" = "SARS-CoV-2 N1", "n1n2" = "SARS-CoV-2 N1N2", "n2" = "SARS-CoV-2 N2")) + 
	scale_fill_manual(name = "Target", values = TARGET_FILLS, labels = c("n1" = "SARS-CoV-2 N1", "n1n2" = "SARS-CoV-2 N1N2", "n2" = "SARS-CoV-2 N2")) + 
	plot_theme() + 
	geom_point(aes(x = day, y = day_sum), color = "#4b1f6f", shape = 1, size = 4, alpha=0.7) + 
	geom_line(aes(x = day, y = day_sum), color = "#4b1f6f")
