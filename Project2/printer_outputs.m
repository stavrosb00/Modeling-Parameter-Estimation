function fig = printer_outputs(t_span, y, y_hat)
    fig = figure;
    plot(t_span, y, 'Linewidth', 1);
    hold on;
    plot(t_span, y_hat, 'Linewidth', 1);
    legend({'$y$', '$\hat{y}$'}, 'Interpreter', 'latex');
    xlabel('$t(sec)$', 'interpreter', 'latex', 'FontWeight', 'bold');
    ylabel('Actual system and simulated model output ');
    
end