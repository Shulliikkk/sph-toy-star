#include "Visualisation.h"

Visualisation::Visualisation() {}

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

void Visualisation::loop(Vec_2d<Particle> result, Vec_2d<double> density, double h) {
    int N_part = result[0].size();
    std::vector<sf::CircleShape> sprites(N_part);
    std::vector<double> rhos;
    rhos.clear();

    for (std::size_t i = 0; i < result.size(); i++) {
        for(std::size_t j = 0; j < result[i].size(); j++) {
            double rho = density[i][j];
            rhos.push_back(rho);
        }
    }

    sf::RenderWindow window (sf::VideoMode (600, 600), "Toy Star Simulation");
    double want_fps = 60;
    sf::Clock loop_timer;

    auto it = std::minmax_element(rhos.begin(), rhos.end(), std::greater<double>());
    auto max_rho = *it.first;
    auto min_rho = *it.second;

    double rhos_range = max_rho - min_rho;
    double rhos_step = rhos_range / 360;

    double rhos_quarter = min_rho + 360 * rhos_step / 4;
    double rhos_half = rhos_quarter + 180 * rhos_step / 2;

    std::sort(rhos.begin(), rhos.end() - 1);

    std::size_t step = 0;
    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
        }

        window.clear();

        for (int i = 0; i < sprites.size(); i++) {
            sprites[i].setRadius(2);
            sprites[i].setOutlineThickness(1.5);

            for (std::size_t j = 0; j < result[step].size(); j++) {
                sprites[i].setPosition(result[step][j].get_position()[0] * 100 + 300, result[step][j].get_position()[1] * 100 + 300);
                double rho = density[step][j];
                if (rho < rhos_quarter)
                {
                    for (int k = 1; k < 360; k++)
                    {
                        if (min_rho + (k - 1) * rhos_step / 4 <= rho && rho < min_rho + k * rhos_step / 4)
                        {
                            sprites[i].setFillColor(hsv(k, 1.f, 0.5));
                            sprites[i].setOutlineColor(hsv(k, 1.f, 1.f));
                        }
                    }
                }

                else if (rho < rhos_half)
                {
                    for (int k = 1; k < 180; k++)
                    {
                        if (rhos_quarter + (k - 1) * rhos_step / 2 <= rho && rho < rhos_quarter + k * rhos_step / 2)
                        {
                            sprites[i].setFillColor(hsv(k, 1.f, 0.5));
                            sprites[i].setOutlineColor(hsv(k, 1.f, 1.f));
                        }
                    }
                }

                else
                {
                    for (int k = 1; k < 180; k++)
                    {
                        if (rhos_half + (k - 1) * rhos_step <= rho && rho < rhos_half + k * rhos_step)
                        {
                            sprites[i].setFillColor(hsv(k, 1.f, 0.5));
                            sprites[i].setOutlineColor(hsv(k, 1.f, 1.f));
                        }
                    }
                }

                window.draw(sprites[i]);
            }
        }

        //шкала
        sf::VertexArray quad1(sf::Quads, 4);
        sf::VertexArray quad2(sf::Quads, 4);
        sf::VertexArray quad3(sf::Quads, 4);

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

       sf::Text text_min;
       sf::Text text_max;
       sf::Text text_half;
       sf::Text text_quarter;

       sf::Font font;
       font.loadFromFile("arial.ttf");

       text_min.setFont(font);
       text_max.setFont(font);
       text_half.setFont(font);
       text_quarter.setFont(font);

       std::string min_str = std::to_string(min_rho);
       std::string max_str = std::to_string(max_rho);
       std::string half_str = std::to_string(rhos_half);
       std::string quarter_str = std::to_string(rhos_quarter);

       text_min.setString(min_str);
       text_max.setString(max_str);
       text_half.setString(half_str);
       text_quarter.setString(quarter_str);

       int size = 12;

       text_min.setCharacterSize(size);
       text_max.setCharacterSize(size);
       text_half.setCharacterSize(size);
       text_quarter.setCharacterSize(size);

       text_min.setFillColor(sf::Color::White);
       text_max.setFillColor(sf::Color::White);
       text_half.setFillColor(sf::Color::White);
       text_quarter.setFillColor(sf::Color::White);

       text_min.setPosition(560.f - 8 * min_str.length(), 442.f - size / 2);
       text_max.setPosition(560.f - 8 * max_str.length(), 152.f - size / 2);
       text_half.setPosition(560.f - 8 * half_str.length(), 297.f - size / 2);
       text_quarter.setPosition(560.f - 8 * quarter_str.length(), 369.5 - size / 2);

       sf::Vertex line1[] =
       {
           sf::Vertex(sf::Vector2f(560.f, 155.f)),
           sf::Vertex(sf::Vector2f(590.f, 155.f))
       };

       sf::Vertex line2[] =
       {
           sf::Vertex(sf::Vector2f(590.f, 155.f)),
           sf::Vertex(sf::Vector2f(590.f, 445.f))
       };

       sf::Vertex line3[] =
       {
           sf::Vertex(sf::Vector2f(590.f, 445.f)),
           sf::Vertex(sf::Vector2f(560.f, 445.f))
       };

       sf::Vertex line4[] =
       {
           sf::Vertex(sf::Vector2f(560.f, 445.f)),
           sf::Vertex(sf::Vector2f(560.f, 155.f))
       };

       sf::Vertex line5[] =
       {
           sf::Vertex(sf::Vector2f(558.f, 300.f)),
           sf::Vertex(sf::Vector2f(590.f, 300.f))
       };

       sf::Vertex line6[] =
       {
           sf::Vertex(sf::Vector2f(558.f, 372.5)),
           sf::Vertex(sf::Vector2f(590.f, 372.5))
       };

       sf::Vertex line7[] =
       {
           sf::Vertex(sf::Vector2f(558.f, 155.f)),
           sf::Vertex(sf::Vector2f(590.f, 155.f))
       };

       sf::Vertex line8[] =
       {
           sf::Vertex(sf::Vector2f(558.f, 445.f)),
           sf::Vertex(sf::Vector2f(590.f, 445.f))
       };

       sf::Vertex line9[] =
       {
           sf::Vertex(sf::Vector2f(558.f, 299.f)),
           sf::Vertex(sf::Vector2f(590.f, 299.f))
       };

       sf::Vertex line10[] =
       {
           sf::Vertex(sf::Vector2f(558.f, 371.5)),
           sf::Vertex(sf::Vector2f(590.f, 371.5))
       };

       sf::Vertex line11[] =
       {
           sf::Vertex(sf::Vector2f(558.f, 154.f)),
           sf::Vertex(sf::Vector2f(590.f, 154.f))
       };

       sf::Vertex line12[] =
       {
           sf::Vertex(sf::Vector2f(558.f, 444.f)),
           sf::Vertex(sf::Vector2f(590.f, 444.f))
       };
       sf::Vertex line13[] =
       {
           sf::Vertex(sf::Vector2f(558.f, 301.f)),
           sf::Vertex(sf::Vector2f(590.f, 301.f))
       };

       sf::Vertex line14[] =
       {
           sf::Vertex(sf::Vector2f(558.f, 373.5)),
           sf::Vertex(sf::Vector2f(590.f, 373.5))
       };

       sf::Vertex line15[] =
       {
           sf::Vertex(sf::Vector2f(558.f, 156.f)),
           sf::Vertex(sf::Vector2f(590.f, 156.f))
       };

       sf::Vertex line16[] =
       {
           sf::Vertex(sf::Vector2f(558.f, 446.f)),
           sf::Vertex(sf::Vector2f(590.f, 446.f))
       };

       double graphic_step = 290 / 360;

       sf::Vertex line17[] =
       {
           sf::Vertex(sf::Vector2f(560.f, 445.f - graphic_step / 4)),
           sf::Vertex(sf::Vector2f(590.f, 445.f - graphic_step / 4))
       };

       sf::Vertex line18[] =
       {
           sf::Vertex(sf::Vector2f(560.f, 445.f - graphic_step / 2)),
           sf::Vertex(sf::Vector2f(590.f, 445.f - graphic_step / 2))
       };
       sf::Vertex line19[] =
       {
           sf::Vertex(sf::Vector2f(560.f, 372.5 - graphic_step / 2)),
           sf::Vertex(sf::Vector2f(590.f, 372.5 - graphic_step / 2))
       };

       sf::Vertex line20[] =
       {
           sf::Vertex(sf::Vector2f(560.f, 372.5 - graphic_step)),
           sf::Vertex(sf::Vector2f(590.f, 372.5 - graphic_step))
       };

       sf::Vertex line21[] =
       {
           sf::Vertex(sf::Vector2f(560.f, 300.f - graphic_step)),
           sf::Vertex(sf::Vector2f(590.f, 300.f - graphic_step))
       };

       sf::Vertex line22[] =
       {
           sf::Vertex(sf::Vector2f(560.f, 300.f - 2 * graphic_step)),
           sf::Vertex(sf::Vector2f(590.f, 300.f - 2 * graphic_step))
       };

       window.draw(text_min);
       window.draw(text_max);
       window.draw(text_half);
       window.draw(text_quarter);

       window.draw(quad1);
       window.draw(quad2);
       window.draw(quad3);
       window.draw(line1, 2, sf::Lines);
       window.draw(line2, 2, sf::Lines);
       window.draw(line3, 2, sf::Lines);
       window.draw(line4, 2, sf::Lines);
       window.draw(line5, 2, sf::Lines);
       window.draw(line6, 2, sf::Lines);
       window.draw(line7, 2, sf::Lines);
       window.draw(line8, 2, sf::Lines);
       window.draw(line9, 2, sf::Lines);
       window.draw(line10, 2, sf::Lines);
       window.draw(line11, 2, sf::Lines);
       window.draw(line12, 2, sf::Lines);
       window.draw(line13, 2, sf::Lines);
       window.draw(line14, 2, sf::Lines);
       window.draw(line15, 2, sf::Lines);
       window.draw(line16, 2, sf::Lines);
       window.draw(line17, 2, sf::Lines);
       window.draw(line18, 2, sf::Lines);
       window.draw(line19, 2, sf::Lines);
       window.draw(line20, 2, sf::Lines);
       window.draw(line21, 2, sf::Lines);
       window.draw(line22, 2, sf::Lines);

       window.display();

        if (step + 1  < result.size())
            step += 1;
        else
            step = 0;

        sf::Int32 frame_duration = loop_timer.getElapsedTime().asMilliseconds();
        sf::Int32 time_to_sleep = int(1000.f/want_fps) - frame_duration;

        if (time_to_sleep > 0)
        {
            sf::sleep(sf::milliseconds(time_to_sleep));
        }

        loop_timer.restart();
    }
}

