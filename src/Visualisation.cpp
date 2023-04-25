#include "Visualisation.h"

sf::Color Visualisation::hsv(int hue, double sat, double val) {
      hue %= 360;
      while (hue < 0) hue += 360;

      if (sat < 0.f) sat = 0.f;
      if (sat > 1.f) sat = 1.f;

      if (val < 0.f) val = 0.f;
      if (val > 1.f) val = 1.f;

      int h = hue / 60;
      double f = double(hue) / 60 - h;
      double p = val * (1.f - sat);
      double q = val * (1.f - sat * f);
      double t = val * (1.f - sat * (1 - f));

      switch(h) {
          default:
          case 0:
          case 6: return sf::Color(val * 255, t * 255, p * 255);
          case 1: return sf::Color(q * 255, val * 255, p * 255);
          case 2: return sf::Color(p * 255, val * 255, t * 255);
          case 3: return sf::Color(p * 255, q * 255, val * 255);
          case 4: return sf::Color(t * 255, p * 255, val * 255);
          case 5: return sf::Color(val * 255, p * 255, q * 255);
      }
}

Visualisation::Visualisation(Vec_2d<Particle> result,
                             Vec_2d<double> density,
                             double h, double system_radius) :
        result(result), density(density), N_part(result[0].size()),
        sprites(N_part),
        h(h), system_radius(system_radius),
        quad1(sf::Quads, 4), quad2(sf::Quads, 4), quad3(sf::Quads, 4) {

    //calc density for scale
    for (auto& step : density) {
        auto bounds = std::minmax_element(step.begin(), step.end());
        rho_steps_bounds.push_back({*bounds.first, *bounds.second});
        /*
        for (std::size_t j = 0; j < N_part; j++) {
            double rho = density[i][j];
            rhos.push_back(rho);
        }
        */
    }

    //std::sort(rhos.begin(), rhos.end());

    //scale creation
    quad1[0].position = sf::Vector2f(560.f, 155.f);
    quad1[1].position = sf::Vector2f(590.f, 155.f);
    quad1[2].position = sf::Vector2f(590.f, 255.f);
    quad1[3].position = sf::Vector2f(560.f, 255.f);

    quad2[0].position = sf::Vector2f(560.f, 255.f);
    quad2[1].position = sf::Vector2f(590.f, 255.f);
    quad2[2].position = sf::Vector2f(590.f, 350.f);
    quad2[3].position = sf::Vector2f(560.f, 350.f);

    quad3[0].position = sf::Vector2f(560.f, 350.f);
    quad3[1].position = sf::Vector2f(590.f, 350.f);
    quad3[2].position = sf::Vector2f(590.f, 445.f);
    quad3[3].position = sf::Vector2f(560.f, 445.f);

    quad1[0].color = sf::Color{255, 0, 255, 255};
    quad1[1].color = sf::Color{255, 0, 255, 255};
    quad1[2].color = sf::Color{0, 255, 255, 255};
    quad1[3].color = sf::Color{0, 255, 255, 255};

    quad2[0].color = sf::Color{0, 255, 255, 255};
    quad2[1].color = sf::Color{0, 255, 255, 255};
    quad2[2].color = sf::Color{255, 255, 0, 255};
    quad2[3].color = sf::Color{255, 255, 0, 255};

    quad3[0].color = sf::Color{255, 255, 0, 255};
    quad3[1].color = sf::Color{255, 255, 0, 255};
    quad3[2].color = sf::Color{255, 0, 0, 255};
    quad3[3].color = sf::Color{255, 0, 0, 255};

    font.loadFromFile("data/arial.ttf");

    text_min.setFont(font);
    text_max.setFont(font);
    text_half.setFont(font);
    text_quarter.setFont(font);
    text_density.setFont(font);

    text_density.setString("Density");

    int size = 14;

    text_min.setCharacterSize(size);
    text_max.setCharacterSize(size);
    text_half.setCharacterSize(size);
    text_quarter.setCharacterSize(size);
    text_density.setCharacterSize(20);

    text_min.setFillColor(sf::Color::White);
    text_max.setFillColor(sf::Color::White);
    text_half.setFillColor(sf::Color::White);
    text_quarter.setFillColor(sf::Color::White);
    text_density.setFillColor(sf::Color::White);

    /*
    text_min.setPosition(555.f - 8 * min_str.length(), 442.f - size / 2);
    text_max.setPosition(555.f - 8 * max_str.length(), 152.f - size / 2);
    text_half.setPosition(555.f - 8 * half_str.length(),
                          297.f - size / 2);
    text_quarter.setPosition(555.f - 8 * quarter_str.length(),
                             369.5 - size / 2);
    text_density.setPosition((570.f + (560.f - 7 * max_str.length())) / 2 - 25.f, 120.f);
    */
    text_min.setPosition(555.f - 8 * 5, 442.f - size / 2);
    text_max.setPosition(555.f - 8 * 5, 152.f - size / 2);
    text_half.setPosition(555.f - 8 * 5, 297.f - size / 2);
    text_quarter.setPosition(555.f - 8 * 5, 369.5 - size / 2);
    text_density.setPosition(520.f, 120.f);

    lines = {
        {sf::Vertex(sf::Vector2f(560.f, 155.f)),
        sf::Vertex(sf::Vector2f(590.f, 155.f))},

        {sf::Vertex(sf::Vector2f(590.f, 155.f)),
        sf::Vertex(sf::Vector2f(590.f, 445.f))},

        {sf::Vertex(sf::Vector2f(590.f, 445.f)),
        sf::Vertex(sf::Vector2f(560.f, 445.f))},

        {sf::Vertex(sf::Vector2f(560.f, 445.f)),
        sf::Vertex(sf::Vector2f(560.f, 155.f))},

        {sf::Vertex(sf::Vector2f(558.f, 300.f)),
        sf::Vertex(sf::Vector2f(590.f, 300.f))},

        {sf::Vertex(sf::Vector2f(558.f, 372.5)),
        sf::Vertex(sf::Vector2f(590.f, 372.5))},

        {sf::Vertex(sf::Vector2f(558.f, 155.f)),
        sf::Vertex(sf::Vector2f(590.f, 155.f))},

        {sf::Vertex(sf::Vector2f(558.f, 445.f)),
        sf::Vertex(sf::Vector2f(590.f, 445.f))},

        {sf::Vertex(sf::Vector2f(558.f, 299.f)),
        sf::Vertex(sf::Vector2f(590.f, 299.f))},

        {sf::Vertex(sf::Vector2f(558.f, 371.5)),
        sf::Vertex(sf::Vector2f(590.f, 371.5))},

        {sf::Vertex(sf::Vector2f(558.f, 444.f)),
        sf::Vertex(sf::Vector2f(590.f, 444.f))},

        {sf::Vertex(sf::Vector2f(558.f, 444.f)),
        sf::Vertex(sf::Vector2f(590.f, 444.f))},

        {sf::Vertex(sf::Vector2f(558.f, 301.f)),
        sf::Vertex(sf::Vector2f(590.f, 301.f))},

        {sf::Vertex(sf::Vector2f(558.f, 373.5)),
        sf::Vertex(sf::Vector2f(590.f, 373.5))},

        {sf::Vertex(sf::Vector2f(558.f, 156.f)),
        sf::Vertex(sf::Vector2f(590.f, 156.f))},

        {sf::Vertex(sf::Vector2f(558.f, 446.f)),
        sf::Vertex(sf::Vector2f(590.f, 446.f))}
    };
}

void Visualisation::round_string(std::string& string) {
    string.erase(std::next(string.begin(), string.find(".") + 2),
                 string.end()
    );
};

void Visualisation::loop() {
    sf::RenderWindow window (sf::VideoMode (600, 600),
                             "Toy Star Simulation"
    );
    double want_fps = 30;
    sf::Clock loop_timer;

    std::size_t step = 0;
    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
        }

        window.clear();
//!!!
        min_rho = rho_steps_bounds[step].first;
        max_rho = rho_steps_bounds[step].second;

        double rhos_range = max_rho - min_rho;
        rhos_step = rhos_range / 360;

        rhos_quarter = min_rho + 360 * rhos_step / 4;
        rhos_half = rhos_quarter + 180 * rhos_step / 2;

        std::string min_str = std::to_string(min_rho);
        std::string max_str = std::to_string(max_rho);
        std::string half_str = std::to_string(rhos_half);
        std::string quarter_str = std::to_string(rhos_quarter);

        round_string(min_str);
        round_string(max_str);
        round_string(half_str);
        round_string(quarter_str);

        text_min.setString(min_str);
        text_max.setString(max_str);
        text_half.setString(half_str);
        text_quarter.setString(quarter_str);

        for (int i = 0; i < N_part; i++) {
            sprites[i].setRadius(2);
            sprites[i].setOutlineThickness(1.5);
            auto pos = result[step][i].get_position();
            auto pos_multiplier = 300 / (system_radius * 3);
            sprites[i].setPosition(pos[0] * pos_multiplier + 300,
                                   pos[1] * pos_multiplier + 300
            );

            double rho = density[step][i];
            auto hue_value = 359 * (rho - min_rho) / (max_rho - min_rho);
            sprites[i].setFillColor(hsv((int)hue_value, 1., 0.5));
            sprites[i].setOutlineColor(hsv((int)hue_value, 1., 1.));

            window.draw(sprites[i]);
        }

        window.draw(text_min);
        window.draw(text_max);
        window.draw(text_half);
        window.draw(text_quarter);
        window.draw(text_density);

        window.draw(quad1);
        window.draw(quad2);
        window.draw(quad3);

        for (auto line : lines) {
            sf::Vertex arr_line[] = {line[0], line[1]};
            window.draw(arr_line, 2, sf::Lines);
        }

        window.display();

        if (step + 1 < result.size()) step++;
        else step = 0;

        sf::Int32 frame_duration = loop_timer.getElapsedTime().asMilliseconds();
        sf::Int32 time_to_sleep = int(1000.f/want_fps) - frame_duration;

        if (time_to_sleep > 0) sf::sleep(sf::milliseconds(time_to_sleep));

        loop_timer.restart();
    }
}

