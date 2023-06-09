function fig = printer_error(t_span, error, is_mse)
    fig = figure;
    plot(t_span, error, 'Linewidth', 1);
    if is_mse
        xlabel('$t(sec)$', 'interpreter', 'latex', 'FontWeight', 'bold');
        ylabel('$\big(y - \hat{y}\big)^2$', 'interpreter', 'latex', ...
                                 'FontWeight', 'bold', 'FontSize', 12);
    else
        xlabel('$t(sec)$', 'interpreter', 'latex', 'FontWeight', 'bold');
        ylabel('$\big|\frac{y-\hat{y}}{y}\big|\%$', 'interpreter', ...
                        'latex', 'FontWeight', 'bold', 'FontSize', 16);
    end
end